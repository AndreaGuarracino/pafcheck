use anyhow::{Context, Result};
use clap::{App, Arg};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

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
                .required_unless("info"),
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
            Arg::with_name("info")
                .long("info")
                .help("Show PAF statistics (mapping lengths, identities, CIGAR presence)")
                .takes_value(false)
                .required(false),
        )
        .get_matches();

    let query_fasta_path = matches.value_of("query_fasta");
    let target_fasta_path = matches.value_of("target_fasta");
    let paf_path = matches.value_of("paf").unwrap();
    let error_mode = matches.value_of("error-mode").unwrap();
    let show_info = matches.is_present("info");

    if show_info {
        if let Err(e) = show_paf_info(paf_path) {
            eprintln!("[pafcheck] Error: {}", e);
            std::process::exit(1);
        }
    } else {
        let query_fasta = query_fasta_path.unwrap();
        let target_fasta = target_fasta_path.unwrap_or(query_fasta);
        if let Err(e) = validate_paf(query_fasta, target_fasta, paf_path, error_mode) {
            eprintln!("[pafcheck] Error: {}", e);
            std::process::exit(1);
        }
    }
}

fn show_paf_info(paf_path: &str) -> Result<()> {
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    let mut total_records = 0;
    let mut records_with_cigar = 0;
    let mut query_lengths = Vec::new();
    let mut target_lengths = Vec::new();
    let mut gap_compressed_identities = Vec::new();
    let mut block_identities = Vec::new();
    
    // Coverage tracking
    let mut query_coverage: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut target_coverage: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut query_sequence_lengths: HashMap<String, usize> = HashMap::new();
    let mut target_sequence_lengths: HashMap<String, usize> = HashMap::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read PAF line")?;
        let record = PafRecord::from_line(&line).context(format!(
            "Failed to parse PAF record at line {}",
            line_number + 1
        ))?;

        total_records += 1;
        
        if !record.cigar.is_empty() {
            records_with_cigar += 1;
        }

        let query_mapping_length = record.query_end - record.query_start;
        let target_mapping_length = record.target_end - record.target_start;
        
        query_lengths.push(query_mapping_length);
        target_lengths.push(target_mapping_length);

        let gap_compressed_identity = if record.block_length > 0 {
            (record.matching_bases as f64 / record.block_length as f64) * 100.0
        } else {
            0.0
        };
        gap_compressed_identities.push(gap_compressed_identity);

        let query_block_identity = if query_mapping_length > 0 {
            (record.matching_bases as f64 / query_mapping_length as f64) * 100.0
        } else {
            0.0
        };
        block_identities.push(query_block_identity);

        // Track coverage
        query_coverage.entry(record.query_name.clone())
            .or_insert_with(Vec::new)
            .push((record.query_start, record.query_end));
        target_coverage.entry(record.target_name.clone())
            .or_insert_with(Vec::new)
            .push((record.target_start, record.target_end));
        
        // Track sequence lengths
        query_sequence_lengths.insert(record.query_name, record.query_length);
        target_sequence_lengths.insert(record.target_name, record.target_length);
    }

    if total_records == 0 {
        println!("No PAF records found");
        return Ok(());
    }

    query_lengths.sort_unstable();
    target_lengths.sort_unstable();
    gap_compressed_identities.sort_by(|a, b| a.partial_cmp(b).unwrap());
    block_identities.sort_by(|a, b| a.partial_cmp(b).unwrap());

    println!("PAF Statistics:");
    println!("===============");
    println!("Total records: {}", total_records);
    println!("Records with CIGAR: {} ({:.1}%)", records_with_cigar, 
             (records_with_cigar as f64 / total_records as f64) * 100.0);
    println!();
    
    println!("Query mapping lengths:");
    print_length_stats(&query_lengths);
    println!();
    
    println!("Target mapping lengths:");
    print_length_stats(&target_lengths);
    println!();
    
    println!("Gap-compressed identity:");
    print_identity_stats(&gap_compressed_identities);
    println!();
    
    println!("Block identity (matches/query_length):");
    print_identity_stats(&block_identities);
    println!();

    // Calculate coverage statistics
    let query_coverage_stats = calculate_coverage_stats(&query_coverage, &query_sequence_lengths);
    let target_coverage_stats = calculate_coverage_stats(&target_coverage, &target_sequence_lengths);

    println!("Query coverage:");
    print_coverage_stats(&query_coverage_stats);
    println!();

    println!("Target coverage:");
    print_coverage_stats(&target_coverage_stats);

    Ok(())
}

fn print_length_stats(lengths: &[usize]) {
    if lengths.is_empty() {
        println!("  No data");
        return;
    }
    
    let total: usize = lengths.iter().sum();
    let mean = total as f64 / lengths.len() as f64;
    let median = if lengths.len() % 2 == 0 {
        (lengths[lengths.len() / 2 - 1] + lengths[lengths.len() / 2]) as f64 / 2.0
    } else {
        lengths[lengths.len() / 2] as f64
    };
    
    println!("  Min: {}", lengths[0]);
    println!("  Max: {}", lengths[lengths.len() - 1]);
    println!("  Mean: {:.1}", mean);
    println!("  Median: {:.1}", median);
}

fn print_identity_stats(identities: &[f64]) {
    if identities.is_empty() {
        println!("  No data");
        return;
    }
    
    let total: f64 = identities.iter().sum();
    let mean = total / identities.len() as f64;
    let median = if identities.len() % 2 == 0 {
        (identities[identities.len() / 2 - 1] + identities[identities.len() / 2]) / 2.0
    } else {
        identities[identities.len() / 2]
    };
    
    println!("  Min: {:.2}%", identities[0]);
    println!("  Max: {:.2}%", identities[identities.len() - 1]);
    println!("  Mean: {:.2}%", mean);
    println!("  Median: {:.2}%", median);
}

#[derive(Debug)]
struct CoverageStats {
    coverage_percentages: Vec<f64>,
    mean: f64,
    variance: f64,
    skewness: f64,
    kurtosis: f64,
    total_sequences: usize,
}

fn calculate_coverage_stats(
    coverage_map: &HashMap<String, Vec<(usize, usize)>>,
    sequence_lengths: &HashMap<String, usize>,
) -> CoverageStats {
    let mut coverage_percentages = Vec::new();
    
    for (seq_name, intervals) in coverage_map {
        if let Some(&seq_length) = sequence_lengths.get(seq_name) {
            let covered_bases = calculate_covered_bases(intervals);
            let coverage_pct = (covered_bases as f64 / seq_length as f64) * 100.0;
            coverage_percentages.push(coverage_pct);
        }
    }
    
    // Add sequences with zero coverage
    for (seq_name, &_seq_length) in sequence_lengths {
        if !coverage_map.contains_key(seq_name) {
            coverage_percentages.push(0.0);
        }
    }
    
    if coverage_percentages.is_empty() {
        return CoverageStats {
            coverage_percentages: Vec::new(),
            mean: 0.0,
            variance: 0.0,
            skewness: 0.0,
            kurtosis: 0.0,
            total_sequences: 0,
        };
    }
    
    coverage_percentages.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let mean = coverage_percentages.iter().sum::<f64>() / coverage_percentages.len() as f64;
    let variance = calculate_variance(&coverage_percentages, mean);
    let skewness = calculate_skewness(&coverage_percentages, mean, variance);
    let kurtosis = calculate_kurtosis(&coverage_percentages, mean, variance);
    
    CoverageStats {
        coverage_percentages,
        mean,
        variance,
        skewness,
        kurtosis,
        total_sequences: sequence_lengths.len(),
    }
}

fn calculate_covered_bases(intervals: &[(usize, usize)]) -> usize {
    if intervals.is_empty() {
        return 0;
    }
    
    let mut sorted_intervals = intervals.to_vec();
    sorted_intervals.sort_by(|a, b| a.0.cmp(&b.0));
    
    let mut covered = 0;
    let mut current_end = 0;
    
    for &(start, end) in &sorted_intervals {
        if start >= current_end {
            covered += end - start;
            current_end = end;
        } else if end > current_end {
            covered += end - current_end;
            current_end = end;
        }
    }
    
    covered
}

fn calculate_variance(values: &[f64], mean: f64) -> f64 {
    if values.len() <= 1 {
        return 0.0;
    }
    
    let sum_squared_diff: f64 = values.iter()
        .map(|x| (x - mean).powi(2))
        .sum();
    
    sum_squared_diff / (values.len() - 1) as f64
}

fn calculate_skewness(values: &[f64], mean: f64, variance: f64) -> f64 {
    if values.len() < 3 || variance == 0.0 {
        return 0.0;
    }
    
    let std_dev = variance.sqrt();
    let sum_cubed_diff: f64 = values.iter()
        .map(|x| ((x - mean) / std_dev).powi(3))
        .sum();
    
    sum_cubed_diff / values.len() as f64
}

fn calculate_kurtosis(values: &[f64], mean: f64, variance: f64) -> f64 {
    if values.len() < 4 || variance == 0.0 {
        return 0.0;
    }
    
    let std_dev = variance.sqrt();
    let sum_fourth_diff: f64 = values.iter()
        .map(|x| ((x - mean) / std_dev).powi(4))
        .sum();
    
    (sum_fourth_diff / values.len() as f64) - 3.0  // Excess kurtosis
}

fn print_coverage_stats(stats: &CoverageStats) {
    if stats.coverage_percentages.is_empty() {
        println!("  No data");
        return;
    }
    
    let median = if stats.coverage_percentages.len() % 2 == 0 {
        let mid = stats.coverage_percentages.len() / 2;
        (stats.coverage_percentages[mid - 1] + stats.coverage_percentages[mid]) / 2.0
    } else {
        stats.coverage_percentages[stats.coverage_percentages.len() / 2]
    };
    
    println!("  Total sequences: {}", stats.total_sequences);
    println!("  Min coverage: {:.2}%", stats.coverage_percentages[0]);
    println!("  Max coverage: {:.2}%", stats.coverage_percentages[stats.coverage_percentages.len() - 1]);
    println!("  Mean coverage: {:.2}%", stats.mean);
    println!("  Median coverage: {:.2}%", median);
    println!("  Std deviation: {:.2}%", stats.variance.sqrt());
    println!("  Skewness: {:.3}", stats.skewness);
    println!("  Kurtosis: {:.3}", stats.kurtosis);
}

fn validate_paf(
    query_fasta: &str,
    target_fasta: &str,
    paf_path: &str,
    error_mode: &str,
) -> Result<()> {
    let mut fasta_reader = MultiFastaReader::new(query_fasta, target_fasta)
        .context("Failed to create FASTA readers")?;
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    let mut total_error_count = 0;
    let mut error_type_counts: HashMap<ErrorType, usize> = HashMap::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read PAF line")?;
        let record = PafRecord::from_line(&line).context(format!(
            "Failed to parse PAF record at line {}",
            line_number + 1
        ))?;

        let mut output = Vec::new();
        if let Err(e) = validate_record(&record, &mut fasta_reader, error_mode, &mut output) {
            if let Some(validation_error) = e.downcast_ref::<ValidationError>() {
                for (error_type, error_info) in &validation_error.errors {
                    let count = error_info.count;
                    *error_type_counts.entry(error_type.clone()).or_insert(0) += count;
                    total_error_count += count;
                    println!(
                        "[pafcheck] Error at line {}: {:?}: {}",
                        line_number + 1,
                        error_type,
                        error_info.first_message
                    );
                    if count > 1 {
                        println!("[pafcheck] {:?}: Total occurrences: {}", error_type, count);
                    }
                }
            } else {
                total_error_count += 1;
                println!("[pafcheck] Error at line {}: {}", line_number + 1, e);
            }
        }
    }

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
