# sum_metrics_parallel.py

import pandas as pd
import numpy as np
import glob
import os
from joblib import Parallel, delayed
from tqdm import tqdm


# === Configuration ===
input_dir = "./21/"  # <-- Update this path with chromosome number
output_file = "Results/Sensitivity_7628_9768_ConfusionMetrics_21.csv"

# Define columns
id_cols = ["CHROM", "POS", "ID","REF", "ALT"]
metric_cols = [
    "TP_HomR", "TP_Het", "TP_HomA",
    "FN_HomR", "FN_Het", "FN_HomA",
    "FP_HomR", "FP_Het", "FP_HomA",
    "TN_HomR", "TN_Het", "TN_HomA"
]

# Auto-detect number of CPUs from SLURM or fallback to all CPUs
n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count()))
print(f"Using {n_jobs} CPU cores")

# Get all CSV files
files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
if not files:
    raise FileNotFoundError(f"No .csv files found in {input_dir}")

# Read IDs from first file only
print(f"Reading ID columns from: {files[0]}")
id_df = pd.read_csv(files[0], usecols=id_cols)

# Function to load metrics from one file
def load_metric_array(file):
    df = pd.read_csv(file, usecols=metric_cols)
    return df.to_numpy(dtype=np.int64)

# Load all files in parallel with progress bar
print(f"Loading and processing {len(files)} files...")
arrays = Parallel(n_jobs=n_jobs)(
    delayed(load_metric_array)(f) for f in tqdm(files, desc="Reading files")
)

# Sum the arrays element-wise
print("Summing arrays...")
sum_array = np.sum(arrays, axis=0)

# Build dataframe with summed columns
sum_df = pd.DataFrame(sum_array, columns=[col + "_Sum62" for col in metric_cols])

# Concatenate ID columns and summed metrics
final_df = pd.concat([id_df, sum_df], axis=1)

# Save to CSV
final_df.to_csv(output_file, index=False)
print(f"✅ Done. Output saved to: {output_file}")
