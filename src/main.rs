use anyhow::{Context, Result};
use clap::{App, Arg};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::thread;
use std::sync::Arc;
use indicatif::{ProgressBar, ProgressStyle};

use pafcheck::fasta_reader::MultiFastaReader;
use pafcheck::paf_parser::PafRecord;
use pafcheck::validator::{validate_record, ErrorType, ValidationError};

fn main() {
    let matches = App::new("PAF Validator")
        .version("1.0")
        .author("Your Name")
        .about("Validates PAF CIGAR strings against FASTA files")
        .arg(
            Arg::with_name("query_fasta")
                .short('q')
                .long("query-fasta")
                .value_name("QUERY_FASTA")
                .help("Path to the bgzip-compressed and tabix-indexed query FASTA file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("target_fasta")
                .short('t')
                .long("target-fasta")
                .value_name("TARGET_FASTA")
                .help("Path to the bgzip-compressed and tabix-indexed target FASTA file")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::with_name("paf")
                .short('p')
                .long("paf")
                .value_name("PAF")
                .help("Path to the PAF file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("error-mode")
                .short('e')
                .long("error-mode")
                .value_name("MODE")
                .help("Error handling mode: omit, report")
                .takes_value(true)
                .required(false)
                .default_value("omit"),
        )
        .arg(
            Arg::with_name("threads")
                .short('j')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads to use for validation")
                .takes_value(true)
                .required(false),
        )
        .get_matches();

    let query_fasta_path = matches.value_of("query_fasta").unwrap();
    let target_fasta_path = matches.value_of("target_fasta").unwrap_or(query_fasta_path);
    let paf_path = matches.value_of("paf").unwrap();
    let error_mode = matches.value_of("error-mode").unwrap();
    
    let num_threads = if let Some(threads_str) = matches.value_of("threads") {
        threads_str.parse::<usize>()
            .context("Invalid number of threads")
            .unwrap_or_else(|e| {
                eprintln!("[pafcheck] Error: {}", e);
                std::process::exit(1);
            })
    } else {
        if std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1) >= 4 { 4 } else { 1 }
    };

    if let Err(e) = validate_paf(query_fasta_path, target_fasta_path, paf_path, error_mode, num_threads) {
        eprintln!("[pafcheck] Error: {}", e);
        std::process::exit(1);
    }
}

struct ThreadResult {
    error_count: usize,
    error_type_counts: HashMap<ErrorType, usize>,
    error_messages: Vec<String>,
}

fn validate_paf(
    query_fasta: &str,
    target_fasta: &str,
    paf_path: &str,
    error_mode: &str,
    num_threads: usize,
) -> Result<()> {
    // Read all PAF lines into memory with line numbers
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);
    
    println!("[pafcheck] Reading PAF file...");
    let lines: Result<Vec<_>> = reader.lines().enumerate().map(|(i, line)| {
        line.context("Failed to read PAF line").map(|l| (i + 1, l))
    }).collect();
    let lines = lines?;
    
    if lines.is_empty() {
        println!("[pafcheck] PAF file is empty. No validation needed.");
        return Ok(());
    }
    
    println!("[pafcheck] Processing {} PAF records using {} threads", lines.len(), num_threads);
    
    // Create progress bar
    let progress_bar = Arc::new(ProgressBar::new(lines.len() as u64));
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("[pafcheck] {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} records ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-")
    );
    
    // Determine chunk size
    let chunk_size = (lines.len() + num_threads - 1) / num_threads;
    
    // Process chunks in parallel using scoped threads
    let results: Result<Vec<ThreadResult>> = thread::scope(|s| {
        let handles: Vec<_> = lines.chunks(chunk_size).enumerate().map(|(chunk_idx, chunk)| {
            let chunk = chunk.to_vec();
            let query_fasta = query_fasta.to_string();
            let target_fasta = target_fasta.to_string();
            let error_mode = error_mode.to_string();
            let progress = Arc::clone(&progress_bar);
            
            s.spawn(move || -> Result<ThreadResult> {
                process_chunk(&query_fasta, &target_fasta, &error_mode, chunk, progress, chunk_idx)
            })
        }).collect();
        
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    
    progress_bar.finish_with_message("Validation complete");
    
    let results = results?;
    
    // Aggregate results
    let mut total_error_count = 0;
    let mut error_type_counts: HashMap<ErrorType, usize> = HashMap::new();
    let mut all_error_messages = Vec::new();
    
    for result in results {
        total_error_count += result.error_count;
        for (error_type, count) in result.error_type_counts {
            *error_type_counts.entry(error_type).or_insert(0) += count;
        }
        all_error_messages.extend(result.error_messages);
    }
    
    // Sort error messages by line number (extract line number from message)
    all_error_messages.sort_by(|a, b| {
        let line_a = extract_line_number(a).unwrap_or(0);
        let line_b = extract_line_number(b).unwrap_or(0);
        line_a.cmp(&line_b)
    });
    
    // Print all error messages in order
    for message in all_error_messages {
        println!("{}", message);
    }
    
    // Print summary
    if total_error_count > 0 {
        println!("[pafcheck] PAF validation completed with errors:");
        for (error_type, count) in error_type_counts.iter() {
            println!("[pafcheck]   - {:?}: {} errors", error_type, count);
        }
        println!("[pafcheck] Total errors: {}", total_error_count);
        anyhow::bail!("PAF validation failed with {} errors", total_error_count);
    } else {
        println!("[pafcheck] PAF validation completed successfully. No errors found.");
        Ok(())
    }
}

fn process_chunk(
    query_fasta: &str,
    target_fasta: &str,
    error_mode: &str,
    chunk: Vec<(usize, String)>,
    progress: Arc<ProgressBar>,
    chunk_idx: usize,
) -> Result<ThreadResult> {
    let mut fasta_reader = MultiFastaReader::new(query_fasta, target_fasta)
        .context("Failed to create FASTA readers")?;
    
    let mut error_count = 0;
    let mut error_type_counts: HashMap<ErrorType, usize> = HashMap::new();
    let mut error_messages = Vec::new();
    
    for (i, (line_number, line)) in chunk.iter().enumerate() {
        // Update progress
        if i % 100 == 0 {
            progress.inc(100.min(chunk.len() - i) as u64);
            if chunk_idx == 0 && i % 1000 == 0 {
                // Only the first thread updates the message to avoid flickering
                progress.set_message(format!("Processing record {}", line_number));
            }
        }
        
        let record = PafRecord::from_line(&line).context(format!(
            "Failed to parse PAF record at line {}",
            line_number
        ))?;

        let mut output = Vec::new();
        if let Err(e) = validate_record(&record, &mut fasta_reader, error_mode, &mut output) {
            if let Some(validation_error) = e.downcast_ref::<ValidationError>() {
                for (error_type, error_info) in &validation_error.errors {
                    let count = error_info.count;
                    *error_type_counts.entry(error_type.clone()).or_insert(0) += count;
                    error_count += count;
                    error_messages.push(format!(
                        "[pafcheck] Error at line {}: {:?}: {}",
                        line_number,
                        error_type,
                        error_info.first_message
                    ));
                    if count > 1 {
                        error_messages.push(format!(
                            "[pafcheck] {:?}: Total occurrences: {}",
                            error_type, count
                        ));
                    }
                }
            } else {
                error_count += 1;
                error_messages.push(format!(
                    "[pafcheck] Error at line {}: {}",
                    line_number, e
                ));
            }
        }
    }
    
    // Update progress for remaining items
    let remaining = chunk.len() % 100;
    if remaining > 0 {
        progress.inc(remaining as u64);
    }
    
    Ok(ThreadResult {
        error_count,
        error_type_counts,
        error_messages,
    })
}

fn extract_line_number(message: &str) -> Option<usize> {
    // Extract line number from error message format "[pafcheck] Error at line X:"
    if let Some(start) = message.find("line ") {
        let start = start + 5; // Skip "line "
        if let Some(end) = message[start..].find(':') {
            return message[start..start + end].parse().ok();
        }
    }
    None
}