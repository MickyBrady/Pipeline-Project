import os
import numpy as np
import pandas as pd

# Define the base directory where kallisto output is stored
base_dir = "/home/2025/mbrady9/PipelineProject_Michaela_Brady"

# Define the log file path
log_file_path = os.path.join(base_dir, "PipelineProject.log")

# Ensure the base directory exists
if not os.path.exists(base_dir):
    os.makedirs(base_dir)

# Ensure the log file has a header
if not os.path.exists(log_file_path) or os.stat(log_file_path).st_size == 0:
    with open(log_file_path, "w") as log_file:
        log_file.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

# Define your samples and conditions
samples_and_conditions = [
    ('Donor1', '2dpi'),
    ('Donor1', '6dpi'),
    ('Donor3', '2dpi'),
    ('Donor3', '6dpi')
]

# Open the log file for appending results
with open(log_file_path, "a") as log_file:
    for sample, condition in samples_and_conditions:
        abundance_file = os.path.join(base_dir, f"Kallisto_Sleuth/{sample}_{condition}", "abundance.tsv")

        if not os.path.exists(abundance_file):
            log_file.write(f"Error: {sample}_{condition} - File not found\n")
            continue

        try:
            df = pd.read_csv(abundance_file, sep="\t")

            if "tpm" not in df.columns:
                log_file.write(f"Error: {sample}_{condition} - TPM column not found\n")
                continue

            tpm_values = df["tpm"].dropna().to_numpy()

            if len(tpm_values) == 0:
                log_file.write(f"Error: {sample}_{condition} - No valid TPM values found\n")
                continue

            min_tpm = np.min(tpm_values)
            med_tpm = np.median(tpm_values)
            mean_tpm = np.mean(tpm_values)
            max_tpm = np.max(tpm_values)

            result_line = f"{sample}\t{condition}\t{min_tpm:.2f}\t{med_tpm:.2f}\t{mean_tpm:.2f}\t{max_tpm:.2f}\n"
            log_file.write(result_line)

        except Exception as e:
            log_file.write(f"Error processing {sample}_{condition}: {e}\n")

print(f"TPM results have been written to {log_file_path}")
