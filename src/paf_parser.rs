use anyhow::{Context, Result};

#[derive(Debug)]
pub struct PafRecord {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub target_name: String,
    pub target_length: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub matching_bases: usize,
    pub block_length: usize,
    pub mapping_quality: usize,
    pub cigar: String,
}

impl PafRecord {
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            anyhow::bail!("PAF line does not have enough fields");
        }

        let mut cigar = String::new();
        for field in &fields[12..] {
            if let Some(stripped) = field.strip_prefix("cg:Z:") {
                cigar = stripped.to_string();
                break;
            }
        }

        Ok(PafRecord {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().context("Failed to parse query length")?,
            query_start: fields[2].parse().context("Failed to parse query start")?,
            query_end: fields[3].parse().context("Failed to parse query end")?,
            strand: fields[4].chars().next().unwrap(),
            target_name: fields[5].to_string(),
            target_length: fields[6].parse().context("Failed to parse target length")?,
            target_start: fields[7].parse().context("Failed to parse target start")?,
            target_end: fields[8].parse().context("Failed to parse target end")?,
            matching_bases: fields[9].parse().context("Failed to parse matching bases")?,
            block_length: fields[10].parse().context("Failed to parse block length")?,
            mapping_quality: fields[11].parse().context("Failed to parse mapping quality")?,
            cigar,
        })
    }
}
