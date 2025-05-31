import os
import re
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Extract TE_AVE and ERROR from last line of files across temperature dirs.")
parser.add_argument("input_filename", help="Name of the input file to read in each subdirectory (e.g., energy.dat)")
parser.add_argument("output_filename", help="Name of the output file to write results (e.g., output.csv)")
args = parser.parse_args()

# Configuration
base_path = "./RESULTS H=0 G"
target_filename = args.input_filename
output_filename = args.output_filename

# Match directories like: out 0.5 h=0 MRT
dir_pattern = re.compile(r"(\d+p\d+)")

data_rows = []

for dirname in os.listdir(base_path):
    full_path = os.path.join(base_path, dirname)
    if not os.path.isdir(full_path):
        continue

    match = dir_pattern.match(dirname)
    if not match:
        continue  # Skip dirs that don't match the pattern

    temperature = float(match.group(1).replace("p", "."))
    file_path = os.path.join(full_path, target_filename)

    if not os.path.isfile(file_path):
        print(f"Skipping missing file: {file_path}")
        continue

    with open(file_path, 'r') as f:
        lines = [line for line in f if line.strip() and not line.strip().startswith("#")]

    if not lines:
        print(f"No data in file: {file_path}")
        continue

    last_line = lines[-1].split()
    try:
        te_ave = float(last_line[-2])
        error = float(last_line[-1])
        data_rows.append((temperature, te_ave, error))
    except (IndexError, ValueError):
        print(f"Could not parse last line in {file_path}: {' '.join(last_line)}")

# Sort by temperature
data_rows.sort()

# Write to output file (CSV style)
with open(output_filename, "w") as out:
    for temp, te_ave, error in data_rows:
        out.write(f"{temp:.3f},{te_ave:.6e},{error:.6e}\n")

print(f"Saved results to {output_filename}")