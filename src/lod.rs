//! LOD (Limit of Detection) calculation and detectability scoring

use crate::{
    bam::process_variant_chunk, DetectabilityResult, LodConfig, Variant, VlodError, VlodResult,
};
use rayon::prelude::*;
use std::path::Path;

/// Chunk variants for parallel processing
pub fn chunkify<T: Clone>(items: Vec<T>, num_chunks: usize) -> Vec<Vec<T>> {
    if items.is_empty() || num_chunks == 0 {
        return vec![items];
    }

    let num_chunks = std::cmp::min(num_chunks, items.len());
    let chunk_size = std::cmp::max(1, items.len() / num_chunks);
    
    let mut chunks = Vec::new();
    let mut start = 0;
    
    for i in 0..num_chunks {
        let end = if i == num_chunks - 1 {
            items.len() // Last chunk gets all remaining items
        } else {
            std::cmp::min(start + chunk_size, items.len())
        };
        
        if start < items.len() {
            chunks.push(items[start..end].to_vec());
            start = end;
        }
    }
    
    chunks
}

/// Calculate detectability scores for a list of variants
pub fn calculate_detectability_scores(
    variants: Vec<Variant>,
    bam_path: &Path,
    config: &LodConfig,
    num_processes: usize,
) -> VlodResult<Vec<DetectabilityResult>> {
    if variants.is_empty() {
        return Ok(Vec::new());
    }

    let num_processes = std::cmp::min(num_processes, variants.len());
    let chunks = chunkify(variants, num_processes);

    // Process chunks in parallel
    let chunk_results: Result<Vec<Vec<_>>, VlodError> = chunks
        .into_par_iter()
        .map(|chunk| process_variant_chunk(&chunk, bam_path, config))
        .collect();

    let chunk_results = chunk_results?;
    
    // Flatten results
    let mut results: Vec<(Variant, f64, u32, u32)> = Vec::new();
    for chunk_result in chunk_results {
        results.extend(chunk_result);
    }

    if results.is_empty() {
        return Ok(Vec::new());
    }

    // Calculate normalization factors (currently unused but kept for potential future use)
    let _max_coverage = results.iter().map(|(_, _, coverage, _)| *coverage).max().unwrap_or(1);
    let _max_variant_reads = results.iter().map(|(_, _, _, reads)| *reads).max().unwrap_or(1);

    // Convert to DetectabilityResult
    let detectability_results: Vec<DetectabilityResult> = results
        .into_iter()
        .map(|(variant, lod, coverage, variant_reads)| {
            let detectability_score = if lod == f64::NEG_INFINITY || coverage <= 1 {
                0.0
            } else {
                lod
            };

            let detectability_condition = if detectability_score >= 2.50 {
                "Detectable".to_string()
            } else {
                "Non-detectable".to_string()
            };

            DetectabilityResult::new(
                variant,
                detectability_score,
                detectability_condition,
                coverage,
                variant_reads,
            )
        })
        .collect();

    Ok(detectability_results)
}

/// Calculate LOD score for a given VAF and configuration
pub fn calculate_lod_score(vaf: f64, config: &LodConfig) -> f64 {
    if vaf <= 0.0 {
        return f64::NEG_INFINITY;
    }

    let lod_value = (config.p_tp * vaf) / ((1.0 - vaf) * config.p_se + vaf * config.p_fp);
    
    if lod_value > 0.0 {
        lod_value.log10()
    } else {
        f64::NEG_INFINITY
    }
}

/// Calculate detectability condition based on score
pub fn calculate_detectability_condition(score: f64) -> String {
    if score >= 2.50 {
        "Detectable".to_string()
    } else {
        "Non-detectable".to_string()
    }
}

/// Validate LOD configuration parameters
pub fn validate_lod_config(config: &LodConfig) -> VlodResult<()> {
    if config.p_tp <= 0.0 || config.p_tp > 1.0 {
        return Err(VlodError::InvalidConfig(
            "p_tp must be between 0 and 1".to_string(),
        ));
    }
    
    if config.p_fp < 0.0 || config.p_fp >= 1.0 {
        return Err(VlodError::InvalidConfig(
            "p_fp must be between 0 and 1".to_string(),
        ));
    }
    
    if config.p_se < 0.0 || config.p_se >= 1.0 {
        return Err(VlodError::InvalidConfig(
            "p_se must be between 0 and 1".to_string(),
        ));
    }

    if config.p_tp <= config.p_fp {
        return Err(VlodError::InvalidConfig(
            "p_tp must be greater than p_fp".to_string(),
        ));
    }

    Ok(())
}

/// Write detectability results to a TSV file
pub fn write_detectability_results(
    results: &[DetectabilityResult],
    output_path: &Path,
) -> VlodResult<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::fs::File;
    use std::io::Write;

    let file = File::create(output_path)?;
    let mut writer: Box<dyn Write> = if output_path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };

    // Write header
    writeln!(
        writer,
        "Chrom\tPos\tRef\tAlt\tDetectability_Score\tDetectability_Condition\tCoverage\tVariant_Reads"
    )?;

    // Write results
    for result in results {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            result.variant.chrom,
            result.variant.pos,
            result.variant.ref_allele,
            result.variant.alt_allele,
            result.detectability_score,
            result.detectability_condition,
            result.coverage,
            result.variant_reads,
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chunkify() {
        let items = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let chunks = chunkify(items, 3);
        
        assert_eq!(chunks.len(), 3);
        assert!(!chunks[0].is_empty());
        assert!(!chunks[1].is_empty());
        assert!(!chunks[2].is_empty());
        
        let total_items: usize = chunks.iter().map(|c| c.len()).sum();
        assert_eq!(total_items, 10);
    }

    #[test]
    fn test_chunkify_empty() {
        let items: Vec<i32> = vec![];
        let chunks = chunkify(items, 3);
        
        assert_eq!(chunks.len(), 1);
        assert!(chunks[0].is_empty());
    }

    #[test]
    fn test_calculate_lod_score() {
        let config = LodConfig::default();
        
        // Test with positive VAF
        let score = calculate_lod_score(0.5, &config);
        assert!(score.is_finite());
        assert!(score > 0.0);
        
        // Test with zero VAF
        let score = calculate_lod_score(0.0, &config);
        assert_eq!(score, f64::NEG_INFINITY);
        
        // Test with negative VAF
        let score = calculate_lod_score(-0.1, &config);
        assert_eq!(score, f64::NEG_INFINITY);
    }

    #[test]
    fn test_calculate_detectability_condition() {
        assert_eq!(calculate_detectability_condition(3.0), "Detectable");
        assert_eq!(calculate_detectability_condition(2.5), "Detectable");
        assert_eq!(calculate_detectability_condition(2.49), "Non-detectable");
        assert_eq!(calculate_detectability_condition(0.0), "Non-detectable");
        assert_eq!(calculate_detectability_condition(-1.0), "Non-detectable");
    }

    #[test]
    fn test_validate_lod_config() {
        let valid_config = LodConfig::default();
        assert!(validate_lod_config(&valid_config).is_ok());
        
        let invalid_config = LodConfig {
            p_tp: 0.0,
            p_fp: 0.001,
            p_se: 0.0001,
        };
        assert!(validate_lod_config(&invalid_config).is_err());
        
        let invalid_config = LodConfig {
            p_tp: 0.5,
            p_fp: 0.6,
            p_se: 0.0001,
        };
        assert!(validate_lod_config(&invalid_config).is_err());
    }
}