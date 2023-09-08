#!/usr/bin/python3

import os
import vcf

# Read files
input_dir = "input_vcf" #input directory containing vcf files
cytoband_file = "../cytoband.txt"  # Path to the cytoband file
output_dir = "out"  # Specify your output directory

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

def process_vcf(vcf_file, cytoband_info, output_dir):
    # Create an output filename based on the input VCF file
    vcf_filename = os.path.basename(vcf_file)
    output_filename = os.path.splitext(vcf_filename)[0] + ".tsv"
    output_file_path = os.path.join(output_dir, output_filename)

    with open(output_file_path, "w") as output_file:
        # Write the header
        output_file.write("ChrNbr\tLength\tPloidy\tCytoband\n")

        # Parse the VCF file
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))

        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS
            end = record.INFO['END']
            ploidy = float(record.INFO['CN'])

            # Find the cytoband information for the current chrom, pos, and end
            cytobands = cytoband_info.get(chrom, [])
            cytoband_str = ""
            for start, stop, band in cytobands:
                if start <= pos <= stop or start <= end <= stop:
                    cytoband_str += band

            # Write to the output file
            output_file.write(f"{chrom}\t{end - pos}\t{ploidy:.1f}\t{chrom}{cytoband_str}\n")

    print(f"Processed: {vcf_file} => Output file saved to: {output_file_path}")

# Read cytoband.txt and store the information in a dictionary
cytoband_info = {}
with open(cytoband_file, "r") as cytoband_f:
    for line in cytoband_f:
        parts = line.strip().split()
        if len(parts) >= 4:
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            band = parts[3]
            cytoband_info[chrom] = cytoband_info.get(chrom, []) + [(start, end, band)]

# Process all VCF files in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".vcf"):
        vcf_file_path = os.path.join(input_dir, filename)
        process_vcf(vcf_file_path, cytoband_info, output_dir)

