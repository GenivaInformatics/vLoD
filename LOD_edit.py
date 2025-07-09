import os
import pysam
import math
import gzip
from concurrent.futures import ProcessPoolExecutor
import argparse


def process_vcf_chunk(args):
    chunk, bam_file_path, bam_index_file_path, p_se, p_fp, p_tp = args
    results = []

    bam_file = pysam.AlignmentFile(
        bam_file_path, "rb", index_filename=bam_index_file_path
    )

    for line in chunk:
        fields = line.strip().split("\t")
        chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
        # vcf_id = f"{chrom}\t{pos}\t{ref}\t{alt}"

        pileup = bam_file.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            stepper="nofilter",
            min_base_quality=0,
            max_depth=1000000,
        )

        alt_alleles = alt.split(",")
        allele_counts = {allele: 0 for allele in alt_alleles}
        ref_count = 0

        for pileup_column in pileup:
            for pileup_read in pileup_column.pileups:
                if not pileup_read.is_refskip:
                    ref_len, alt_len = (
                        len(ref),
                        max([len(allele) for allele in alt_alleles]),
                    )
                    if ref_len == alt_len:
                        if ref_len == 1:
                            # SNV
                            if not pileup_read.is_del:
                                base = pileup_read.alignment.query_sequence[
                                    pileup_read.query_position
                                ]
                                if base == ref:
                                    ref_count += 1
                                elif base in allele_counts:
                                    allele_counts[base] += 1
                        else:
                            # MNV
                            if not pileup_read.is_del:
                                start = pileup_read.query_position
                                end = start + ref_len
                                read_seq = pileup_read.alignment.query_sequence[
                                    start:end
                                ]
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
            lod = math.log10(lod_value) if lod_value > 0 else float("-inf")

            # vcf_id_alt = f"{chrom}\t{pos}\t{ref}\t{alt}"
            coverage = total_count
            num_variant_reads = alt_count

            results.append((chrom, pos, ref, alt, lod, coverage, num_variant_reads))

    return results


def chunkify(lst, num_chunks):
    # Fix: Ensure chunk_size is at least 1 and num_chunks doesn't exceed list length
    num_chunks = min(num_chunks, len(lst))
    if num_chunks <= 0:
        return [lst]  # Return the whole list as a single chunk if num_chunks is invalid

    chunk_size = max(1, len(lst) // num_chunks)
    return [lst[i : i + chunk_size] for i in range(0, len(lst), chunk_size)]


def is_gzipped(file_path):
    with open(file_path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detectability script")
    parser.add_argument("--input-vcf", help="Path to the input VCF file", required=True)
    parser.add_argument("--input-bam", help="Path to the input BAM file", required=True)
    parser.add_argument(
        "--input-bam-index",
        help="Path to the index file for the BAM file",
        required=True,
    )
    parser.add_argument("--output", help="Path to the output TSV file", required=True)
    parser.add_argument(
        "--TP", type=float, help="Probability of true positive result", default=0.999
    )
    parser.add_argument(
        "--FP", type=float, help="Probability of false positive result", default=0.001
    )
    parser.add_argument(
        "--SE", type=float, help="Probability of sequencing error", default=0.0001
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        help="Number of processes to use",
        default=os.cpu_count(),
    )

    args = parser.parse_args()

    vcf_file_path = args.input_vcf
    bam_file_path = args.input_bam
    bam_index_file_path = args.input_bam_index
    output_file_path = args.output
    num_processes = args.num_processes

    p_tp = args.TP
    p_fp = args.FP
    p_se = args.SE

    vcf_lines = []
    with (
        gzip.open(vcf_file_path, "rt")
        if is_gzipped(vcf_file_path)
        else open(vcf_file_path)
    ) as vcf_file:
        for line in vcf_file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                try:
                    int(fields[1])
                    vcf_lines.append(line)
                except ValueError:
                    continue

    # Adjust num_processes if there are fewer lines than requested processes
    num_processes = min(num_processes, len(vcf_lines)) if len(vcf_lines) > 0 else 1
    num_chunks = num_processes

    # Use the fixed chunkify function
    chunks = chunkify(vcf_lines, num_chunks)
    args = [
        (chunk, bam_file_path, bam_index_file_path, p_se, p_fp, p_tp)
        for chunk in chunks
    ]

    # Handle the case when there are no vcf_lines
    if not vcf_lines:
        results = []
    else:
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            results = list(executor.map(process_vcf_chunk, args))

    results = [result for chunk_results in results for result in chunk_results]

    # Check if results list is empty
    if not results:
        print("No variants found in the input VCF file.")
        with gzip.open(output_file_path, "wt") as output_file:
            output_file.write(
                "Chrom\tPos\tRef\tAlt\tDetectability_Score\tDetectability_Condition\tCoverage\tVariant_Reads\n"
            )
        exit(0)

    # max_lod = max([result[4] for result in results])
    max_coverage = max([result[5] for result in results])
    max_num_variant_reads = max([result[6] for result in results])

    w_lod, w_coverage, w_num_variant_reads = 1, 1, 1

    detectability_results = []
    for chrom, pos, ref, alt, lod, coverage, num_variant_reads in results:
        if lod == float("-inf") or coverage <= 1:
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

        detectability_results.append(
            (chrom, pos, ref, alt, detectability_score, detectability_condition)
        )

    output_data = []
    for (
        chrom,
        pos,
        ref,
        alt,
        detectability_score,
        detectability_condition,
    ) in detectability_results:
        # Get coverage and num_variant_reads for this variant
        for result in results:
            if (
                result[0] == chrom
                and result[1] == pos
                and result[2] == ref
                and result[3] == alt
            ):
                coverage, num_variant_reads = result[5], result[6]
                break

        output_data.append(
            [
                chrom,
                pos,
                ref,
                alt,
                detectability_score,
                detectability_condition,
                coverage,
                num_variant_reads,
            ]
        )

    with gzip.open(output_file_path, "wt") as output_file:
        output_file.write(
            "Chrom\tPos\tRef\tAlt\tDetectability_Score\tDetectability_Condition\tCoverage\tVariant_Reads\n"
        )
        for line in output_data:
            output_file.write("\t".join([str(x) for x in line]) + "\n")
