import os
import sys
import pysam
import math
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import argparse

def process_vcf_chunk(args):
    chunk, bam_file_path, bam_index_file_path, p_se, p_fp, p_tp = args
    results = []

    bam_file = pysam.AlignmentFile(bam_file_path, "rb", index_filename=bam_index_file_path)

    for line in chunk:
        fields = line.strip().split("\t")
        chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
        vcf_id = f"{chrom}_{pos}_{ref}_{alt}"

        pileup = bam_file.pileup(chrom, pos - 1, pos, truncate=True, stepper="nofilter", min_base_quality=0, max_depth=1000000)

        alt_alleles = alt.split(',')
        allele_counts = {allele: 0 for allele in alt_alleles}
        ref_count = 0

        for pileup_column in pileup:
            for pileup_read in pileup_column.pileups:
                if not pileup_read.is_refskip:
                    ref_len, alt_len = len(ref), max([len(allele) for allele in alt_alleles])
                    if ref_len == alt_len:
                        if ref_len == 1:
                            # SNV
                            if not pileup_read.is_del:
                                base = pileup_read.alignment.query_sequence[pileup_read.query_position]
                                if base == ref:
                                    ref_count += 1
                                elif base in allele_counts:
                                    allele_counts[base] += 1
                        else:
                            # MNV
                            if not pileup_read.is_del:
                                start = pileup_read.query_position
                                end = start + ref_len
                                read_seq = pileup_read.alignment.query_sequence[start:end]
                                if read_seq == ref:
                                    ref_count += 1
                                elif read_seq in allele_counts:
                                    allele_counts[read_seq] += 1
                    else:
                        # Indel
                        indel = pileup_read.indel
                        for alt_allele in alt_alleles:
                            if indel == len(alt_allele) - len(ref):
                                allele_counts[alt_allele] += 1
                            elif indel == 0:
                                ref_count += 1

        # Aggregate results for multiple alternative alleles
        total_count = ref_count + sum(allele_counts.values())
        for alt, alt_count in allele_counts.items():
            vaf = alt_count / total_count if total_count > 0 else 0
            lod_value = (p_tp * vaf) / ((1 - vaf) * p_se + vaf * p_fp)
            lod = math.log10(lod_value) if lod_value > 0 else float('-inf')

            vcf_id_alt = f"{chrom}_{pos}_{ref}_{alt}"
            coverage = total_count
            num_variant_reads = alt_count

            results.append((vcf_id_alt, lod, coverage, num_variant_reads))

    return results

def chunkify(lst, num_chunks):
    chunk_size = len(lst) // num_chunks
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Detectability script")
    parser.add_argument("--input-vcf", help="Path to the input VCF file", required=True)
    parser.add_argument("--input-bam", help="Path to the input BAM file", required=True)
    parser.add_argument("--input-bam-index", help="Path to the index file for the BAM file", required=True)
    parser.add_argument("--output", help="Path to the output TSV file", required=True)
    parser.add_argument("--TP", type=float, help="Probability of true positive result", default=0.999)
    parser.add_argument("--FP", type=float, help="Probability of false positive result", default=0.001)
    parser.add_argument("--SE", type=float, help="Probability of sequencing error", default=0.0001)

    args = parser.parse_args()

    vcf_file_path = args.input_vcf
    bam_file_path = args.input_bam
    bam_index_file_path = args.input_bam_index
    output_file_path = args.output

    p_tp = args.TP
    p_fp = args.FP
    p_se = args.SE

    vcf_lines = []
    with open(vcf_file_path, "r") as vcf_file:
        for line in vcf_file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                try:
                    int(fields[1])
                    vcf_lines.append(line)
                except ValueError:
                    continue

    num_chunks = os.cpu_count() or 1
    chunks = chunkify(vcf_lines, num_chunks)
    args = [(chunk, bam_file_path, bam_index_file_path, p_se, p_fp, p_tp) for chunk in chunks]

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_vcf_chunk, args))

    results = [result for chunk_results in results for result in chunk_results]
    max_lod = max([result[1] for result in results])
    max_coverage = max([result[2] for result in results])
    max_num_variant_reads = max([result[3] for result in results])

    w_lod, w_coverage, w_num_variant_reads = 1, 1, 1

    detectability_results = []
    for vcf_id, lod, coverage, num_variant_reads in results:
        if lod == float('-inf') or coverage <= 1:
            detectability_condition = "Non-detectable"
            detectability_score = 0
        else:
            detectability_score = lod
            normalized_coverage = coverage / max_coverage
            normalized_num_variant_reads = num_variant_reads / max_num_variant_reads

        if detectability_score >= 2.50:
            detectability_condition = "Detectable"
        elif detectability_score > 0.00 and detectability_score < 2.50:
            detectability_condition = "Non-detectable"
        else:
            detectability_condition = "Non-detectable"

        detectability_results.append((vcf_id, detectability_score, detectability_condition))

    output_data = []
    for vcf_id, detectability_score, detectability_condition, coverage, num_variant_reads in [(vcf_id, detectability_score, detectability_condition) + (result[2], result[3]) for (vcf_id, detectability_score, detectability_condition), result in zip(detectability_results, results)]:
        if detectability_condition == "Non-detectable":
            output_data.append([vcf_id, detectability_score, detectability_condition, coverage, num_variant_reads])
        else:
            output_data.append([vcf_id, detectability_score, detectability_condition, coverage, num_variant_reads])

    output_df = pd.DataFrame(output_data, columns=["VCF_ID", "Detectability_Score", "Detectability_Condition", "Coverage", "Variant_Reads"])
    output_df.to_csv(output_file_path, index=False, sep="\t")
