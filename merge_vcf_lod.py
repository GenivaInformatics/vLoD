import pandas as pd
import sys
import gzip

def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


def merge_detectability_into_vcf(vcf_path, detectability_path, output_path):
    # Read the detectability status file
    tsv_df = pd.read_csv(detectability_path, sep='\t')
    detectability_data_dict = {}
    for index, row in tsv_df.iterrows():
        detectability_data_dict[(row['Chrom'], row['Pos'], row['Ref'], row['Alt'])] = {
            'condition': 'Yes' if row['Detectability_Condition'] == 'Detectable' else 'No',
            'score': row['Detectability_Score']
        }

    # print( detectability_data_dict)

    info_added = False
    # Process the VCF file
    updated_vcf_lines = []
    with (gzip.open(vcf_path, 'rt') if is_gzipped(vcf_path) else open(vcf_path)) as vcf_file:
        for line in vcf_file:
            if line.startswith("#C"):
                header = line.strip().split()
                info_idx = header.index("INFO")
            if line.startswith("##INFO"):
                updated_vcf_lines.append(line)
                # if "TLOD" in line:
                if not info_added:
                    updated_vcf_lines.append('##INFO=<ID=DET,Number=1,Type=String,Description="Detectability status (Yes if detectable, No if non-detectable)">\n')
                    updated_vcf_lines.append('##INFO=<ID=DETS,Number=1,Type=Float,Description="Detectability Score">\n')
                    info_added = True
                continue
            if line.startswith("##") or line.startswith("#"):
                updated_vcf_lines.append(line)
                continue
            columns = line.strip().split("\t")
            # vcf_id = f"{columns[0]}_{columns[1]}_{columns[3]}_{columns[4]}"
            vcf_id = (columns[0], int(columns[1]), columns[3], columns[4])
            # print(vcf_id)
            if vcf_id in detectability_data_dict:
                detectability_status = detectability_data_dict[vcf_id]['condition']
                detectability_score = detectability_data_dict[vcf_id]['score']
                columns[info_idx] = columns[info_idx] + f";DET={detectability_status};DETS={detectability_score}"
                updated_vcf_lines.append("\t".join(columns) + "\n")
            else:
                updated_vcf_lines.append(line)

    # Save the updated VCF
    with open(output_path, "w") as out_vcf:
        out_vcf.writelines(updated_vcf_lines)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_detectability.py <path_to_vcf> <path_to_detectability_status> <output_path>")
        sys.exit(1)

    vcf_path = sys.argv[1]
    detectability_path = sys.argv[2]
    output_path = sys.argv[3]
    merge_detectability_into_vcf(vcf_path, detectability_path, output_path)