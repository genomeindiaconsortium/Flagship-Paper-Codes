import os
import pandas as pd
import logging
import subprocess

# Setup logging
logging.basicConfig(filename='vep_processing.log', level='DEBUG', format='%(asctime)s %(levelname)s:%(message)s')

def modify_vcf_header(input_vcf):
    """Modify the VCF header to handle any issues with the '#CHROM' line."""
    try:
        with open(input_vcf, 'r') as vcf_file:
            vcf_lines = vcf_file.readlines()

        modified_vcf_lines = ['#' + line if line.startswith('#CHROM') else line for line in vcf_lines]

        with open(input_vcf, 'w') as vcf_file:
            vcf_file.writelines(modified_vcf_lines)

        logging.info(f'Successfully modified VCF header for {input_vcf}')
    except Exception as e:
        logging.error(f'Error modifying VCF header for {input_vcf}: {e}')
        raise

def modify_vcf_header_output(input_vcf):
    """Modify the VCF header to replace '##CHROM' with '#CHROM'."""
    try:
        with open(input_vcf, 'r') as vcf_file:
            vcf_lines = vcf_file.readlines()

        # Replace '##CHROM' with '#CHROM' in the header
        modified_vcf_lines = [line.replace('##CHROM', '#CHROM') if '##CHROM' in line else line for line in vcf_lines]

        with open(input_vcf, 'w') as vcf_file:
            vcf_file.writelines(modified_vcf_lines)

        logging.info(f"Successfully modified VCF header for {input_vcf}")
    except Exception as e:
        logging.error(f"Error modifying VCF header for {input_vcf}: {e}")
        raise








def run_vep(input_vcf, vep_cache, output_dir):
    """Run the VEP tool on the input VCF."""
    try:
        base_name = os.path.splitext(os.path.basename(input_vcf))[0]
        output_file = os.path.join(output_dir, f"{base_name}_output.vcf")
        command = (
             f"vep -i {input_vcf} --cache --dir_cache {vep_cache} --vcf --fork 32 --offline --force_overwrite --flag_pick --everything -o {output_file}  "
             f"--plugin LoF,loftee_path:/data/sauma/.vep/Plugins/loftee-1.0.4_GRCh38/,human_ancestor_fa:/data/sauma/.vep/Plugins/loftee-1.0.4_GRCh38/test/human_ancestor.fa.gz,conservation_file:/data/sauma/.vep/Plugins/loftee-1.0.4_GRCh38/loftee.sql,gerp_bigwig:/data/sauma/.vep/Plugins/loftee-1.0.4_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw  --dir_plugins /data/sauma/.vep/Plugins/loftee-1.0.4_GRCh38/ ")

        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f'VEP command failed for {input_vcf}: {result.stderr}')
        else:
            logging.info(f'Successfully ran VEP for {input_vcf}: {result.stdout}')
        return output_file
    except Exception as e:
        logging.error(f'Error running VEP for {input_vcf}: {e}')
        raise


def generate_csv_from_vcf(input_vcf, output_csv):
    """Generate a CSV file from the VCF file with the desired '~' delimiter and handle multiple annotations."""
    try:
        base_name = os.path.splitext(os.path.basename(input_vcf))[0]
        temp_header = f"{base_name}_temp_header.txt"
        temp_data = f"{base_name}_temp_data.txt"

        # Extract the header format
        header_command = (
            f"bcftools view -h {input_vcf} | grep '##INFO=<ID=CSQ' | "
            f"sed -E 's/.*Format: //' | sed 's/\">//' | tr '|' '~' > {temp_header}"
        )
        os.system(header_command)
        logging.info(f"Generated CSV header for {base_name}.")

        # Read and modify the header to include CHROM~POS
        with open(temp_header, "r") as header_file:
            header = header_file.read().strip()
        header = f"CHROM~POS~{header}"

        # Extract data with replaced delimiters ('|' to '~')
        data_command = (
            f"bcftools query -f '%CHROM|%POS|%INFO/CSQ\\n' {input_vcf} | sed 's/|/~/g' > {temp_data}"
        )
        os.system(data_command)
        logging.info(f"Extracted data with '~' as the delimiter for {base_name}.")

        # Process data to handle multiple annotations
        with open(output_csv, "w") as output_file:
            output_file.write(header + "\n")  # Write the header
            with open(temp_data, "r") as data_file:
                for line in data_file:
                    split_line = line.strip().split("~", 2)
                    if len(split_line) < 3:
                        logging.warning(f"Skipping malformed line: {line.strip()}")
                        continue
                    chrom, pos, annotations = split_line
                    chrom_pos = f"{chrom}~{pos}"
                    for annotation in annotations.split(","):
                        output_file.write(f"{chrom_pos}~{annotation}\n")
        logging.info(f"Final CSV saved to {output_csv} for {base_name}")

        # Clean up unique temp files
        os.remove(temp_header)
        os.remove(temp_data)

    except Exception as e:
        logging.error(f"Error generating CSV from VCF {input_vcf}: {e}")
        raise

def process_vcf(input_vcf, vep_cache, output_dir):
    """Process a single VCF file from header modification to VEP output parsing."""
    try:
        modify_vcf_header(input_vcf)
        output_file = run_vep(input_vcf, vep_cache, output_dir)
        output_csv = os.path.join(output_dir, os.path.splitext(os.path.basename(input_vcf))[0] + "_processed.csv")
        modify_vcf_header_output(output_file)
        generate_csv_from_vcf(output_file, output_csv)
    except Exception as e:
        logging.error(f'Error processing VCF {input_vcf}: {e}')

if __name__ == "__main__":
    vep_cache = "/data/sauma/.vep"
    output_dir = os.getcwd()

    # List of input VCF files
    input_vcf_files = [
    "9768Samples_MAC2_MAC3_Final_chr1_v3-1.vcf.nosamples.vcf"
            ]

    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    
    for input_vcf in input_vcf_files:
        process_vcf(input_vcf, vep_cache, output_dir)
