# for CRISPRDetect
import glob
import pandas as pd
import os
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Align.AlignInfo import SummaryInfo
# Define the file pattern
file_pattern = ""

# Create an empty DataFrame to store the combined data
combined_CRISPRDetect_df = pd.DataFrame(columns=["crispr_id", "contig_id", "crispr_start", "crispr_end", "crispr_len", "direct_repeat", "n_repeats", "repeat_len", "spacer_len_avg"])

for filepath in glob.glob(file_pattern):
    with open(filepath, "r") as f:
        print("Processing file:", filepath)
        crispr_num = 0
        crispr_pos_ls = []
        table = False
        ln_count = 0
        for line in f:
            line = line.strip()
            
            fields = [word.strip() for word in line.split(" ") if word.strip()]
            tab_fields = line.split("\t")
            if line.startswith(">"):
                seqname = tab_fields[0].replace(">", "")
                crispr_num += 1
                crispr_g = str(crispr_num)

            if line.startswith("="):
                ln_count += 1
                continue

            if ln_count == 1:
                crispr_pos = int(fields[0])
                crispr_pos_ls.append(crispr_pos)

            if line.startswith("="):
                ln_count += 1
                continue

            if line and line[0].isdigit() and ln_count == 2:
                line_strip = line.strip()
                crispr_pos_ls = crispr_pos_ls.pop()
                n_repeats = tab_fields[0]
                repeat_len = tab_fields[1].strip()
                spacer_len = tab_fields[3]
                repeat_seq = tab_fields[4]
                crispr_start = min([crispr_pos_ls])
                crispr_end = max([crispr_pos_ls]) + int(repeat_len) - 1
                crispr_len = int(crispr_end) - int(crispr_start) + 1
                print(n_repeats,crispr_start, crispr_end)
                row = {
                        "crispr_id": seqname + "_c" + crispr_g,
                        "contig_id": seqname,
                        "crispr_start": crispr_start,
                        "crispr_end": crispr_end,
                        "crispr_len": crispr_len,
                        "direct_repeat": repeat_seq,
                        "n_repeats": n_repeats,
                        "repeat_len": repeat_len,
                        "spacer_len_avg": int(spacer_len)
                    }
            #     print(row)
                    # row = {"crispr_id": crispr_id, "contig_id": seqname, "crispr_start": crispr_start, "crispr_end": crispr_end, "crispr_len": crispr_len, "direct_repeat": repeat_seq, "n_repeats": n_repeats, "repeat_len": repeat_len, "spacer_len_avg": spacer_len}
                combined_CRISPRDetect_df = combined_CRISPRDetect_df.append(row, ignore_index=True)

# Write the combined data to a single output file
output_filename = ""  # Specify the desired output file path
combined_CRISPRDetect_df.to_csv(output_filename, sep=',', index=False)
print(f"Reformatted data written to {output_filename}.")
