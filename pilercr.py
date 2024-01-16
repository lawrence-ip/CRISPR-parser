# For PilerCR
import glob
import pandas as pd
import os

# Define the file pattern
file_pattern = ""

# Create an empty DataFrame to store the combined data
combined_df = pd.DataFrame(columns=["crispr_id", "contig_id", "crispr_start", "crispr_end", "crispr_len", "direct_repeat", "n_repeats", "repeat_len", "spacer_len_avg"])

# Iterate over matching files
for filepath in glob.glob(file_pattern):
    seqname = os.path.splitext(os.path.basename(filepath))[0]
    process_lines = False
    table = False
    proceed = False
    with_report = False
    
    with open(filepath, "r") as f:
        content = f.read()
        if "DETAIL REPORT" in content:
            print("Processing file:", filepath)
            with_report = True
        else:
            print("No CRISPRS found in file:", filepath) 
            
    with open(filepath, "r") as f:
        if "DETAIL REPORT" in f.read():
            print("Processing file:", filepath)
            with_report = True
            f.seek(0)  # Reset the file pointer to the beginning of the file
            
            crispr_start_ls = []
            smallest_crispr_pos = int('9999999999')
            largest_crispr_pos = int('0')
            
            for line in f:
                print(line)
                fields = [word.strip() for word in line.split(" ") if word.strip()]
                
                if line.startswith("Array"):
                    process_lines = True
                    print("Read Array")
                    array_num = line.split()[1]
                    crispr_id = seqname + "_c" + array_num
                    crispr_start_ls = []  # Reset the list for each array
                    smallest_crispr_pos = int('9999999999')  # Reset the smallest position
                    largest_crispr_pos = int('0')  # Reset the largest position
                
                if line.startswith("="):
                    table = True
                
                if process_lines and table and line.startswith("    ") and not line.startswith("       Pos") and not line.startswith("       "):
                    crispr_pos = int(fields[0])
                    crispr_start_ls.append(crispr_pos)
                    if min(crispr_start_ls) < smallest_crispr_pos:
                        smallest_crispr_pos = min(crispr_start_ls)
                    if max(crispr_start_ls) > largest_crispr_pos:
                        largest_crispr_pos = max(crispr_start_ls)
                    
                if line.startswith("="):
                    table = False
                if process_lines and line.startswith(f"    {largest_crispr_pos}"):
                    crispr_end = int(fields[0]) + int(fields[1]) - 1
                if process_lines and not table and line.startswith("    ") and not line.startswith("       Pos") and not fields[3].isdigit() and len(fields[3]) > 12:
                    crispr_len = str(int(crispr_end) - int(smallest_crispr_pos) + 1)
                    n_repeats = fields[0]
                    repeat_len = fields[1]
                    spacer_len = fields[2]
                    repeat_seq = fields[3]

                if line.startswith("SUMMARY"):
                    process_lines = False 
                
                row = {"crispr_id": crispr_id, "contig_id": seqname, "crispr_start": smallest_crispr_pos, "crispr_end": crispr_end, "crispr_len": crispr_len, "direct_repeat": repeat_seq, "n_repeats": n_repeats, "repeat_len": repeat_len, "spacer_len_avg": spacer_len}
                combined_df = combined_df.append(row, ignore_index=True)
# Write the combined data to a single output file
output_filename = ""  # Specify the desired output file path
combined_df.to_csv(output_filename, sep=',', index=False)
print(f"Reformatted data written to {output_filename}.")
