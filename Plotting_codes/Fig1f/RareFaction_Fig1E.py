import os
import random
from collections import defaultdict
import pandas as pd
from mpi4py import MPI


# ============================================================
#  User settings
# ============================================================
input_folder = "Input"
output_folder = "Output"

# Create output folder if it does not already exist
os.makedirs(output_folder, exist_ok=True)

# File containing sample-to-cluster mapping
'''
The file looks like this:

 Sample     Clusters
 S1         AA_T
 S2         AA_T
 S7800      TB_T 

'''
cluster_file = "sample_clusters.tsv"

chr_num = ["22"] 
# Example full chromosome list:
# chr_num = ["1", "2", "3", "4", "5", "6", "7", "8", "9","10", "11", "12", "13", "14", "15", "16", "17","18", "19", "20", "21"]


# Number of rarefaction permutations per cluster
N_ITERATIONS = 100

# For reproducible sample order generation
RANDOM_SEED = 42

# Maximum number of ranks/samples to retain per permutation.
# If None, use all samples available in the cluster.
MAX_RANKS = None

# ============================================================
#  Helper: construct chromosome file path
# ============================================================

def chrom_file_path(chrnum: str) -> str:
    return os.path.join(input_folder, f"chr{chrnum}_summaryData_novel.csv")


# ============================================================
#  Load sample-to-cluster mapping
'''

    Expected input file columns:
        Sample      sample ID
        Clusters    cluster label

    Returns:
        mapping: dictionary of the form
                 {sample_id: cluster_name}
'''

# ============================================================

def load_sample_cluster_map(path: str) -> dict:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cluster mapping file not found: {path}")

    df = pd.read_csv(path, sep="\t", dtype=str)

    required_cols = {"Sample", "Clusters"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(
            f"Cluster file missing columns: {missing}. "
            f"Columns found: {df.columns.tolist()}"
        )

    mapping = {}

    for _, row in df.iterrows():
        sample_id = str(row["Sample"]).strip()
        cluster = str(row["Clusters"]).strip()

        if sample_id and cluster and cluster != "NA":
            mapping[sample_id] = cluster

    print(f"[CLUSTERS] Loaded {len(mapping)} sample→cluster mappings from {path}")

    return mapping

# Load the sample-to-cluster map once
SAMPLE_TO_CLUSTER = load_sample_cluster_map(cluster_file)


# =========================================================
# BUILD FIXED SAMPLE ORDERS PER CLUSTER
'''
    To generate fixed random sample orders for each cluster.

    For each cluster:
        - collects all samples belonging to that cluster
        - generates n_iterations random orderings
        - truncates each ordering to max_ranks

    Returns:
        fixed_orders[cluster_name][iteration] = [sample1, sample2, ...]

    The same sample orders ensures comparability and aggregation of novel counts across chromosomes
'''

# =========================================================

def build_fixed_cluster_orders(
    sample_to_cluster: dict,
    n_iterations: int = 100,
    seed: int = 42,
    max_ranks=None
):
    cluster_to_all_samples = defaultdict(list)
 
    # Group samples by cluster
    for sample_id, cluster in sample_to_cluster.items():
        cluster_to_all_samples[cluster].append(sample_id)

    fixed_orders = {}

    for cluster_name, sample_ids in cluster_to_all_samples.items():

        # Sorting gives stable input before random shuffling
        sample_ids = sorted(sample_ids)

        n_samples = len(sample_ids)
        n_ranks = n_samples if max_ranks is None else min(max_ranks, n_samples)

        # Base random generator
        base_rng = random.Random(seed)

        cluster_orders = {}

        for it in range(1, n_iterations + 1):

            # Create an iteration-specific random generator
            # using a seed drawn from the base generator.
            rng = random.Random(base_rng.randint(1, 10**9))

            order = sample_ids[:]
            rng.shuffle(order)

             # Keep only the desired number of ranks/samples
            order = order[:n_ranks]

            cluster_orders[it] = order

        fixed_orders[cluster_name] = cluster_orders
    return fixed_orders

# Generate fixed rarefaction sample orders for all clusters
FIXED_CLUSTER_ORDERS = build_fixed_cluster_orders(
    sample_to_cluster=SAMPLE_TO_CLUSTER,
    n_iterations=N_ITERATIONS,
    seed=RANDOM_SEED,
    max_ranks=MAX_RANKS
)

# ============================================================
# Parse sample-variant relationships from one chromosome file
'''
    Parse a chromosome-level variant file and construct:

        sample_to_variants:
            sample_id -> set of variants carried by that sample

        variant_to_samples:
            variant -> set of samples carrying that variant

    Expected columns:
        'Concordance' : variant identifier (chr:pos:ref:alt)
        'DR_Tribe', 'IE_NonTribe', 'DR_NonTribe', 'AA_Tribe', 'TB_NonTribe', 'TB_Tribe', 'IE_Tribe':
            group-specific columns containing sample information

    The function expects each group-column entry to contain colon-separated
    fields, where the 8th field contains semicolon-separated sample entries.
    Refer to file named:  chr22_toy_input.csv 
'''
# ============================================================

def parse_sample_variants_from_chr_file(file_path, label=""):
    if not os.path.exists(file_path):
        print(f"[{label}] File not found: {file_path}")
        return None, None

    try:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"[{label}] Error reading file: {e}")
        return None, None

    if "Concordance" not in df.columns:
        print(f"[{label}] Missing 'Concordance' column. Columns: {df.columns.tolist()}")
        return None, None

    expected_group_cols = [
        "DR_Tribe", "IE_NonTribe", "DR_NonTribe",
        "AA_Tribe", "TB_NonTribe", "TB_Tribe", "IE_Tribe"
    ]

    group_cols = [c for c in expected_group_cols if c in df.columns]

    if not group_cols:
        print(f"[{label}] No expected group columns found.")
        print(f"[{label}] Available columns: {df.columns.tolist()}")
        return None, None

    print(f"[{label}] Using group columns: {group_cols}")

    sample_to_variants = defaultdict(set)
    variant_to_samples = defaultdict(set)
 
    # Loop over variants
    for _, row in df.iterrows():
        variant = row["Concordance"]
        if pd.isna(variant) or str(variant).strip() == "":
            continue
        variant = str(variant).strip()

        # Search all group columns for samples carrying this variant
        for gcol in group_cols:
            col_val = row.get(gcol, None)
            if col_val is None or pd.isna(col_val) or col_val == "NA":
                continue
            
            # Expected structure:
            # field1:field2:...:field8
            # where field8 contains sample entries
            parts = str(col_val).split(":", 7)
            if len(parts) < 8:
                continue
            samples_str = parts[7]
            if not samples_str:
                continue
            
            # Multiple sample entries are separated by semicolons   
            for entry in samples_str.split(";"):
                entry = entry.strip()
                if not entry:
                    continue
                # Each entry may contain sample_id|additional_info
                # Keep only the sample ID before "|"
                sample_id = entry.split("|", 1)[0].strip()
                if not sample_id:
                    continue
                sample_to_variants[sample_id].add(variant)
                variant_to_samples[variant].add(sample_id)
    return sample_to_variants, variant_to_samples


# ============================================================
# Compute rarefaction curves for one cluster
'''
    Inputs:
        sample_to_vars:
            dictionary mapping each sample to the set of variants it carries

        fixed_orders:
            dictionary mapping permutation number to ordered sample list

    For each permutation:
        - start with an empty variant set
        - add samples one-by-one according to the fixed order
        - after each sample, record the cumulative number of unique variants seen

    Returns:
        count_df:
            rows = permutations
            columns = rank_1, rank_2, ...
            values = cumulative number of unique variants discovered

        sample_df:
            rows = permutations
            columns = rank_1, rank_2, ...
            values = sample IDs used at each rank
'''

# ============================================================

def compute_permutation_outputs(
    sample_to_vars,
    fixed_orders,
    label=""
):
    if not fixed_orders:
        print(f"[{label}] No fixed orders available.")
        return None, None
    count_rows = []
    sample_rows = []

    total_permutations = len(fixed_orders)

    for it, order in fixed_orders.items():
        if it == 1 or it % 10 == 0 or it == total_permutations:
            print(f"[{label}]  {it}/{total_permutations}")
        seen = set()
        cumulative_curve = []

         # Add samples sequentially according to the fixed rarefaction order
        for sid in order:

            seen.update(sample_to_vars.get(sid, set()))

            # Store cumulative number of unique variants discovered so far
            cumulative_curve.append(len(seen))

        count_row = {"Permutation": it}
        sample_row = {"Permutation": it}

        for idx, sid in enumerate(order):
            col = f"rank_{idx + 1}"

            count_row[col] = cumulative_curve[idx]
            sample_row[col] = sid

        count_rows.append(count_row)
        sample_rows.append(sample_row)

    count_df = pd.DataFrame(count_rows)
    sample_df = pd.DataFrame(sample_rows)

    return count_df, sample_df

# ============================================================
#  Process one chromosome
'''
Steps:
        1. Read chromosome-level variant file.
        2. Parse sample-variant mappings.
        3. For each cluster:
              a) subset to samples from that cluster
              b) compute rarefaction curves
              c) save cumulative count matrix
              d) save sample-order matrix
'''
# ============================================================

def process_chromosome(chrnum: str):

    label = f"CHR {chrnum}"
    file_path = chrom_file_path(chrnum)

    print(f"[{label}] Reading: {file_path}")

    sample_to_variants, variant_to_samples = parse_sample_variants_from_chr_file(
        file_path,
        label=label
    )

    if sample_to_variants is None:
        print(f"[{label}] No sample-variant mapping returned.")
        return

    print(
        f"[{label}] Parsed "
        f"{len(sample_to_variants)} samples and "
        f"{len(variant_to_samples)} variants."
    )


    # Process each cluster separately   
    for cluster_name, fixed_orders in FIXED_CLUSTER_ORDERS.items():

        cluster_label = f"{label} | {cluster_name}"

        # Collect all samples that appear in any fixed order for this cluster
        cluster_samples = set()

        for order in fixed_orders.values():
            cluster_samples.update(order)

        # Keep only cluster samples that have variants on this chromosome
        sub_map = {
            sid: sample_to_variants[sid]
            for sid in cluster_samples
            if sid in sample_to_variants and sample_to_variants[sid]
        }

        if not sub_map:
            print(
                f"[{cluster_label}] No samples with variants on this chromosome. "
            )

        count_df, sample_df = compute_permutation_outputs(
            sample_to_vars=sub_map,
            fixed_orders=fixed_orders,
            label=cluster_label
        )

        if count_df is None or sample_df is None:
            print(f"[{cluster_label}] No output generated.")
            continue

        # Save cumulative direct-count matrix for one chromosome
        out_count = os.path.join(
            output_folder,
            f"chr{chrnum}_cluster_{cluster_name}_direct_counts.tsv"
        )

        count_df.to_csv(out_count, sep="\t", index=False)

        print(f"[{cluster_label}] Saved direct-count matrix: {out_count}")

         # Save sample names used at each rank for each permutation
        out_samples = os.path.join(
            output_folder,
            f"chr{chrnum}_cluster_{cluster_name}_sample_names.tsv"
        )

        sample_df.to_csv(out_samples, sep="\t", index=False)

        print(f"[{cluster_label}] Saved sample-name matrix: {out_samples}")

# ============================================================
#  Parallelize across chromosomes
# ============================================================

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank < len(chr_num):
    assigned_chr = chr_num[rank]
    print(f"[MPI rank {rank}/{size}] → processing chromosome {assigned_chr}")
    process_chromosome(assigned_chr)
else:
    print(f"[MPI rank {rank}/{size}] → idle (no chromosome assigned)")