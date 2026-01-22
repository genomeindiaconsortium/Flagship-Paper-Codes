import pandas as pd
import glob
import os
import gc
import sys

# ================= CONFIGURATION =================
CLINVAR_FILE = 'variant_summary_all.txt'
CHUNK_SIZE = 100000
OUTDIR = "maf_intermediate"

os.makedirs(OUTDIR, exist_ok=True)
# =================================================

def parse_position_from_snp_id(snp_id):
    try:
        return int(str(snp_id).split(':')[1])
    except Exception:
        return -1

def main():

    # ---------------------------------------------------------
    # 1. Load and Filter ClinVar
    # ---------------------------------------------------------
    print(f"Loading ClinVar: {CLINVAR_FILE}")

    df = pd.read_csv(CLINVAR_FILE, sep='\t', low_memory=False)

    filters = [
        "criteria provided, multiple submitters, no conflicts",
        "reviewed by expert panel",
        "practice guideline"
    ]

    df = df[df["ReviewStatus"].isin(filters)]
    df = df[df["Assembly"] == "GRCh38"].copy()
    
    df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", regex=False)
    df["PositionVCF"] = pd.to_numeric(df["PositionVCF"], errors="coerce")

    df["LookupKey"] = list(zip(
        df["Chromosome"],
        df["PositionVCF"],
        df["ReferenceAlleleVCF"],
        df["AlternateAlleleVCF"]
    ))

    df.set_index("LookupKey", inplace=True)
    valid_keys = set(df.index)

    print(f"ClinVar variants tracked: {len(valid_keys)}")

    # ---------------------------------------------------------
    # 2. Process each FRQ file safely
    # ---------------------------------------------------------
    frq_files = sorted(glob.glob("*.frq"))

    for f_path in frq_files:
        pop = os.path.basename(f_path).replace(".frq", "").replace("_2col", "")
        out_file = f"{OUTDIR}/maf_{pop}.tsv"

        print(f"\nProcessing {pop}")

        # Write header immediately
        with open(out_file, "w") as fh:
            fh.write("Chromosome\tPosition\tRef\tAlt\tMAF\n")

        try:
            chunk_iter = pd.read_csv(
                f_path,
                sep=r"\s+",
                chunksize=CHUNK_SIZE,
                dtype=str
            )
        except Exception as e:
            print(f"❌ Could not read {f_path}: {e}")
            continue

        total_matches = 0

        for chunk in chunk_iter:
            if not {"CHR", "SNP", "A1", "A2", "MAF"}.issubset(chunk.columns):
                continue

            chunk["ParsedPos"] = chunk["SNP"].apply(parse_position_from_snp_id)
            chunk = chunk[chunk["ParsedPos"] != -1]

            rows_to_write = []

            for chr_v, pos_v, a1, a2, maf in zip(
                chunk["CHR"],
                chunk["ParsedPos"],
                chunk["A1"],
                chunk["A2"],
                chunk["MAF"]
            ):
                k1 = (chr_v, pos_v, a2, a1)
                k2 = (chr_v, pos_v, a1, a2)

                if k1 in valid_keys:
                    rows_to_write.append((chr_v, pos_v, a2, a1, maf))
                elif k2 in valid_keys:
                    rows_to_write.append((chr_v, pos_v, a1, a2, maf))

            if rows_to_write:
                out_df = pd.DataFrame(
                    rows_to_write,
                    columns=["Chromosome", "Position", "Ref", "Alt", "MAF"]
                )
                out_df.to_csv(
                    out_file,
                    sep="\t",
                    index=False,
                    header=False,
                    mode="a"
                )
                total_matches += len(rows_to_write)

        print(f"✔ {pop}: {total_matches} matches written to {out_file}")
        gc.collect()

    print("\nAll populations processed safely.")

if __name__ == "__main__":
    main()
