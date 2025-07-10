//! VCF integration functionality for merging detectability results

use crate::{
    vcf::is_gzipped,
    DetectabilityResult, VlodError, VlodResult,
};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// Read detectability results from a TSV file
pub fn read_detectability_results<P: AsRef<Path>>(
    path: P,
) -> VlodResult<HashMap<(String, u32, String, String), (String, f64)>> {
    let file = File::open(&path)
        .map_err(|_| VlodError::FileNotFound(path.as_ref().to_string_lossy().to_string()))?;

    let reader: Box<dyn BufRead> = if is_gzipped(&path)? {
        let gz_decoder = MultiGzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);

    let mut detectability_data = HashMap::new();

    for result in csv_reader.records() {
        let record = result?;
        
        if record.len() < 6 {
            continue;
        }

        let chrom = record[0].to_string();
        let pos = record[1].parse::<u32>()
            .map_err(|_| VlodError::InvalidVariant(format!("Invalid position: {}", &record[1])))?;
        let ref_allele = record[2].to_string();
        let alt_allele = record[3].to_string();
        let detectability_score = record[4].parse::<f64>()
            .map_err(|_| VlodError::InvalidVariant(format!("Invalid score: {}", &record[4])))?;
        let detectability_condition = record[5].to_string();

        let condition = if detectability_condition == "Detectable" {
            "Yes".to_string()
        } else {
            "No".to_string()
        };

        detectability_data.insert(
            (chrom, pos, ref_allele, alt_allele),
            (condition, detectability_score),
        );
    }

    Ok(detectability_data)
}

/// Merge detectability results into a VCF file
pub fn merge_detectability_into_vcf<P: AsRef<Path>>(
    vcf_path: P,
    detectability_path: P,
    output_path: P,
) -> VlodResult<()> {
    let detectability_data = read_detectability_results(detectability_path)?;

    let file = File::open(&vcf_path)
        .map_err(|_| VlodError::FileNotFound(vcf_path.as_ref().to_string_lossy().to_string()))?;

    let reader: Box<dyn BufRead> = if is_gzipped(&vcf_path)? {
        let gz_decoder = MultiGzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut output_file = File::create(output_path)?;
    let mut info_added = false;
    let mut info_column_index = None;

    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with("#CHROM") {
            // Find the INFO column index
            let header: Vec<&str> = line.split('\t').collect();
            info_column_index = header.iter().position(|&col| col == "INFO");
            writeln!(output_file, "{}", line)?;
            continue;
        }

        if line.starts_with("##INFO") {
            writeln!(output_file, "{}", line)?;
            if !info_added {
                writeln!(
                    output_file,
                    "##INFO=<ID=DET,Number=1,Type=String,Description=\"Detectability status (Yes if detectable, No if non-detectable)\">"
                )?;
                writeln!(
                    output_file,
                    "##INFO=<ID=DETS,Number=1,Type=Float,Description=\"Detectability Score\">"
                )?;
                info_added = true;
            }
            continue;
        }

        if line.starts_with("##") || line.starts_with("#") {
            writeln!(output_file, "{}", line)?;
            continue;
        }

        // Process data lines
        let mut columns: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        
        if columns.len() < 8 {
            writeln!(output_file, "{}", line)?;
            continue;
        }

        let chrom = columns[0].clone();
        let pos = columns[1].parse::<u32>().unwrap_or(0);
        let ref_allele = columns[3].clone();
        let alt_allele = columns[4].clone();

        let vcf_id = (chrom, pos, ref_allele, alt_allele);

        if let Some((condition, score)) = detectability_data.get(&vcf_id) {
            let info_idx = info_column_index.unwrap_or(7);
            
            if info_idx < columns.len() {
                let new_info = format!("{};DET={};DETS={}", columns[info_idx], condition, score);
                columns[info_idx] = new_info;
            }
        }

        writeln!(output_file, "{}", columns.join("\t"))?;
    }

    Ok(())
}

/// Create detectability results from a vector of DetectabilityResult
pub fn create_detectability_map(
    results: &[DetectabilityResult],
) -> HashMap<(String, u32, String, String), (String, f64)> {
    let mut map = HashMap::new();
    
    for result in results {
        let key = (
            result.variant.chrom.clone(),
            result.variant.pos,
            result.variant.ref_allele.clone(),
            result.variant.alt_allele.clone(),
        );
        
        let condition = if result.detectability_condition == "Detectable" {
            "Yes".to_string()
        } else {
            "No".to_string()
        };
        
        map.insert(key, (condition, result.detectability_score));
    }
    
    map
}

/// Merge detectability results directly into VCF without intermediate file
pub fn merge_detectability_results_into_vcf<P: AsRef<Path>>(
    vcf_path: P,
    results: &[DetectabilityResult],
    output_path: P,
) -> VlodResult<()> {
    let detectability_data = create_detectability_map(results);

    let file = File::open(&vcf_path)
        .map_err(|_| VlodError::FileNotFound(vcf_path.as_ref().to_string_lossy().to_string()))?;

    let reader: Box<dyn BufRead> = if is_gzipped(&vcf_path)? {
        let gz_decoder = MultiGzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut output_file = File::create(output_path)?;
    let mut info_added = false;
    let mut info_column_index = None;

    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with("#CHROM") {
            // Find the INFO column index
            let header: Vec<&str> = line.split('\t').collect();
            info_column_index = header.iter().position(|&col| col == "INFO");
            writeln!(output_file, "{}", line)?;
            continue;
        }

        if line.starts_with("##INFO") {
            writeln!(output_file, "{}", line)?;
            if !info_added {
                writeln!(
                    output_file,
                    "##INFO=<ID=DET,Number=1,Type=String,Description=\"Detectability status (Yes if detectable, No if non-detectable)\">"
                )?;
                writeln!(
                    output_file,
                    "##INFO=<ID=DETS,Number=1,Type=Float,Description=\"Detectability Score\">"
                )?;
                info_added = true;
            }
            continue;
        }

        if line.starts_with("##") || line.starts_with("#") {
            writeln!(output_file, "{}", line)?;
            continue;
        }

        // Process data lines
        let mut columns: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        
        if columns.len() < 8 {
            writeln!(output_file, "{}", line)?;
            continue;
        }

        let chrom = columns[0].clone();
        let pos = columns[1].parse::<u32>().unwrap_or(0);
        let ref_allele = columns[3].clone();
        let alt_allele = columns[4].clone();

        let vcf_id = (chrom, pos, ref_allele, alt_allele);

        if let Some((condition, score)) = detectability_data.get(&vcf_id) {
            let info_idx = info_column_index.unwrap_or(7);
            
            if info_idx < columns.len() {
                let new_info = format!("{};DET={};DETS={}", columns[info_idx], condition, score);
                columns[info_idx] = new_info;
            }
        }

        writeln!(output_file, "{}", columns.join("\t"))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Variant;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_detectability_results() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "Chrom\tPos\tRef\tAlt\tDetectability_Score\tDetectability_Condition\tCoverage\tVariant_Reads").unwrap();
        writeln!(temp_file, "chr1\t100\tA\tT\t3.5\tDetectable\t30\t15").unwrap();
        writeln!(temp_file, "chr2\t200\tG\tC\t1.2\tNon-detectable\t20\t5").unwrap();
        
        let results = read_detectability_results(temp_file.path()).unwrap();
        
        assert_eq!(results.len(), 2);
        assert_eq!(results.get(&("chr1".to_string(), 100, "A".to_string(), "T".to_string())), Some(&("Yes".to_string(), 3.5)));
        assert_eq!(results.get(&("chr2".to_string(), 200, "G".to_string(), "C".to_string())), Some(&("No".to_string(), 1.2)));
    }

    #[test]
    fn test_create_detectability_map() {
        let variant = Variant::new("chr1".to_string(), 100, "A".to_string(), "T".to_string());
        let result = DetectabilityResult::new(
            variant,
            3.5,
            "Detectable".to_string(),
            30,
            15,
        );
        
        let map = create_detectability_map(&[result]);
        
        assert_eq!(map.len(), 1);
        assert_eq!(map.get(&("chr1".to_string(), 100, "A".to_string(), "T".to_string())), Some(&("Yes".to_string(), 3.5)));
    }

    #[test]
    fn test_merge_detectability_into_vcf() {
        // Create test detectability file
        let mut detectability_file = NamedTempFile::new().unwrap();
        writeln!(detectability_file, "Chrom\tPos\tRef\tAlt\tDetectability_Score\tDetectability_Condition\tCoverage\tVariant_Reads").unwrap();
        writeln!(detectability_file, "chr1\t100\tA\tT\t3.5\tDetectable\t30\t15").unwrap();
        
        // Create test VCF file
        let mut vcf_file = NamedTempFile::new().unwrap();
        writeln!(vcf_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf_file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
        writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(vcf_file, "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30").unwrap();
        
        let output_file = NamedTempFile::new().unwrap();
        
        merge_detectability_into_vcf(
            vcf_file.path(),
            detectability_file.path(),
            output_file.path(),
        ).unwrap();
        
        // Read the output and verify
        let output_content = std::fs::read_to_string(output_file.path()).unwrap();
        assert!(output_content.contains("DET=Yes"));
        assert!(output_content.contains("DETS=3.5"));
        assert!(output_content.contains("##INFO=<ID=DET,Number=1,Type=String"));
        assert!(output_content.contains("##INFO=<ID=DETS,Number=1,Type=Float"));
    }
}