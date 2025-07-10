//! Utility functions for file handling and common operations

use crate::{VlodError, VlodResult};
use std::fs::File;
use std::io::Read;
use std::path::Path;

/// Check if a file is gzip compressed
pub fn is_gzipped<P: AsRef<Path>>(path: P) -> VlodResult<bool> {
    let mut file = File::open(path)?;
    let mut buffer = [0; 2];
    
    match file.read_exact(&mut buffer) {
        Ok(()) => Ok(buffer == [0x1f, 0x8b]),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(VlodError::Io(e)),
    }
}

/// Get the number of CPU cores, with a fallback default
pub fn get_num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
}

/// Validate file paths and check if they exist
pub fn validate_file_exists<P: AsRef<Path>>(path: P) -> VlodResult<()> {
    if !path.as_ref().exists() {
        return Err(VlodError::FileNotFound(
            path.as_ref().to_string_lossy().to_string(),
        ));
    }
    Ok(())
}

/// Validate that a file is readable
pub fn validate_file_readable<P: AsRef<Path>>(path: P) -> VlodResult<()> {
    validate_file_exists(&path)?;
    
    File::open(&path)
        .map_err(|_| VlodError::FileNotFound(path.as_ref().to_string_lossy().to_string()))?;
    
    Ok(())
}

/// Check if a path has a specific extension
pub fn has_extension<P: AsRef<Path>>(path: P, extension: &str) -> bool {
    path.as_ref()
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.eq_ignore_ascii_case(extension))
        .unwrap_or(false)
}

/// Get file extension as a string
pub fn get_extension<P: AsRef<Path>>(path: P) -> Option<String> {
    path.as_ref()
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_lowercase())
}

/// Format a file size in bytes to a human-readable string
pub fn format_file_size(size: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = size as f64;
    let mut unit_index = 0;
    
    while size >= 1024.0 && unit_index < UNITS.len() - 1 {
        size /= 1024.0;
        unit_index += 1;
    }
    
    format!("{:.2} {}", size, UNITS[unit_index])
}

/// Create parent directories if they don't exist
pub fn ensure_parent_dirs<P: AsRef<Path>>(path: P) -> VlodResult<()> {
    if let Some(parent) = path.as_ref().parent() {
        std::fs::create_dir_all(parent)?;
    }
    Ok(())
}

/// Log progress information
pub fn log_progress(current: usize, total: usize, message: &str) {
    if total > 0 {
        let percentage = (current as f64 / total as f64) * 100.0;
        log::info!("{}: {} / {} ({:.1}%)", message, current, total, percentage);
    }
}

/// Chunking utility for splitting work across threads
pub fn chunk_work<T: Clone>(items: Vec<T>, num_chunks: usize) -> Vec<Vec<T>> {
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

/// Timer utility for measuring execution time
pub struct Timer {
    start: std::time::Instant,
    name: String,
}

impl Timer {
    pub fn new(name: &str) -> Self {
        log::info!("Starting timer: {}", name);
        Timer {
            start: std::time::Instant::now(),
            name: name.to_string(),
        }
    }
    
    pub fn elapsed(&self) -> std::time::Duration {
        self.start.elapsed()
    }
    
    pub fn log_elapsed(&self) {
        let duration = self.elapsed();
        log::info!("Timer '{}' elapsed: {:.2?}", self.name, duration);
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        self.log_elapsed();
    }
}

/// Memory usage reporting utility
pub fn log_memory_usage(context: &str) {
    #[cfg(unix)]
    {
        use std::fs;
        if let Ok(status) = fs::read_to_string("/proc/self/status") {
            for line in status.lines() {
                if line.starts_with("VmRSS:") {
                    if let Some(memory_str) = line.split_whitespace().nth(1) {
                        if let Ok(memory_kb) = memory_str.parse::<u64>() {
                            let memory_mb = memory_kb / 1024;
                            log::info!("Memory usage ({}): {} MB", context, memory_mb);
                        }
                    }
                    break;
                }
            }
        }
    }
    
    #[cfg(not(unix))]
    {
        log::debug!("Memory usage logging not supported on this platform ({})", context);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_is_gzipped() {
        // Test with a regular file
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "test content").unwrap();
        assert_eq!(is_gzipped(temp_file.path()).unwrap(), false);
        
        // Test with gzipped content
        let mut temp_file = NamedTempFile::new().unwrap();
        temp_file.write_all(&[0x1f, 0x8b]).unwrap();
        assert_eq!(is_gzipped(temp_file.path()).unwrap(), true);
    }

    #[test]
    fn test_get_num_cpus() {
        let num_cpus = get_num_cpus();
        assert!(num_cpus >= 1);
    }

    #[test]
    fn test_validate_file_exists() {
        let temp_file = NamedTempFile::new().unwrap();
        assert!(validate_file_exists(temp_file.path()).is_ok());
        
        assert!(validate_file_exists("/nonexistent/file").is_err());
    }

    #[test]
    fn test_has_extension() {
        assert!(has_extension("test.vcf", "vcf"));
        assert!(has_extension("test.VCF", "vcf"));
        assert!(!has_extension("test.txt", "vcf"));
        assert!(!has_extension("test", "vcf"));
    }

    #[test]
    fn test_get_extension() {
        assert_eq!(get_extension("test.vcf"), Some("vcf".to_string()));
        assert_eq!(get_extension("test.VCF"), Some("vcf".to_string()));
        assert_eq!(get_extension("test"), None);
    }

    #[test]
    fn test_format_file_size() {
        assert_eq!(format_file_size(512), "512.00 B");
        assert_eq!(format_file_size(1024), "1.00 KB");
        assert_eq!(format_file_size(1536), "1.50 KB");
        assert_eq!(format_file_size(1024 * 1024), "1.00 MB");
    }

    #[test]
    fn test_chunk_work() {
        let items = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let chunks = chunk_work(items, 3);
        
        assert!(chunks.len() <= 3);
        assert!(!chunks.is_empty());
        
        let total_items: usize = chunks.iter().map(|c| c.len()).sum();
        assert_eq!(total_items, 10);
    }

    #[test]
    fn test_chunk_work_empty() {
        let items: Vec<i32> = vec![];
        let chunks = chunk_work(items, 3);
        
        assert_eq!(chunks.len(), 1);
        assert!(chunks[0].is_empty());
    }

    #[test]
    fn test_timer() {
        let timer = Timer::new("test");
        std::thread::sleep(std::time::Duration::from_millis(1));
        assert!(timer.elapsed().as_millis() >= 1);
    }
}