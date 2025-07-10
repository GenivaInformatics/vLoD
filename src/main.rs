fn main() {
    println!("vlod-rs - Variant Limit of Detection Tool");
    println!();
    println!("🔬 RECOMMENDED: Use the combined tool for most workflows:");
    println!("  vlod          - Complete analysis: VCF + BAM → Annotated VCF (one step)");
    println!();
    println!("📋 Advanced tools for specialized workflows:");
    println!("  lod_edit      - Analyze variant detectability (VCF + BAM → TSV)");
    println!("  merge_vcf_lod - Merge detectability results into VCF (VCF + TSV → VCF)");
    println!();
    println!("📖 For help with each tool:");
    println!("  cargo run --help                        # This combined tool");
    println!("  cargo run --bin lod_edit -- --help      # Analysis only");
    println!("  cargo run --bin merge_vcf_lod -- --help  # Merging only");
    println!();
    println!("🚀 Quick start example:");
    println!("  cargo run -- --input-vcf variants.vcf --input-bam alignments.bam --output annotated.vcf");
    println!();
    println!("💡 The combined tool eliminates intermediate files and is faster than the two-step process.");
}
