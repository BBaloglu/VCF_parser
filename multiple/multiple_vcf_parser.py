#!/usr/bin/python3

import os
import vcf
import glob

# Define the input VCF files, cytoband file, output directory, and output file prefix
#input_vcf_files = ["input1.vcf", "input2.vcf"]  # Add your input VCF file paths here
input_vcf_files = glob.glob('*.vcf')
cytoband_file = "../cytoband.txt"  # Path to the cytoband file
output_directory = "out"  # Specify your output directory

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

def get_cytoband(chrom, pos, end, cytobands):
    for start, stop, band in cytobands:
        if int(start) <= pos <= int(stop) or int(start) <= end <= int(stop):
            return f"{chrom}{band}"
    return "N/A"

# Process each input VCF file
for input_vcf in input_vcf_files:
    # Create a VCF reader
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))

    # Read cytoband information into a list of tuples
    cytobands = []
    with open(cytoband_file, 'r') as cytoband_f:
        for line in cytoband_f:
            parts = line.strip().split('\t')
            cytobands.append((parts[1], parts[2], parts[3]))

    # Create the output file
    output_file_path = os.path.join(output_directory, os.path.splitext(os.path.basename(input_vcf))[0] + ".tsv") 
    
    with open(output_file_path, 'w') as output_file:
        output_file.write("ChrNbr\tLength\tPloidy\tCytoband\n")
        
        # Iterate through VCF records
        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS
            end = record.INFO['END']
            
            # Extract ploidy from the QUAL field
            ploidy = record.QUAL
            
            # Get the cytoband information
            cytoband = get_cytoband(chrom, pos, end, cytobands)
            
            # Write to the output file
            output_file.write(f"{chrom}\t{end - pos}\t{ploidy}\t{cytoband}\n")

print("Processing complete.")

