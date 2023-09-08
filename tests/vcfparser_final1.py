#!/usr/bin/python3
import vcf

def parse_cytoband(cytoband_file):
    cytoband_info = {}
    with open(cytoband_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chromosome = parts[0]
            arm_band = parts[3]
            cytoband_info[chromosome] = arm_band
    return cytoband_info

def extract_cnv_and_generate_output(vcf_file, cytoband_info, output_file):
    with open(vcf_file, 'r') as vcf_reader, open(output_file, 'w') as output:
        # Write the header for the output file
        output.write("ChrNbr\tLength\tPloidy\tCytoband\n")

        for line in vcf_reader:
            if line.startswith('#'):
                continue  # Skip header lines
            parts = line.strip().split('\t')
            chromosome = parts[0]
            end_position = int(parts[7].split('END=')[1].split(';')[0])
            length = end_position - int(parts[1])
            ploidy = float(parts[9].split(':')[2])
            arm_band = cytoband_info.get(chromosome, 'Unknown')

            # Write the information to the output file
            output.write(f"{chromosome}\t{length}\t{ploidy:.2f}\t{arm_band}\n")

# Example usage
cytoband_file = '../cytoband.txt'
vcf_file = 'data/CN_Segments_mosaic_Seq46.vcf'
output_file = 'cytoband_arm_Seq46.tsv'

cytoband_info = parse_cytoband(cytoband_file)
extract_cnv_and_generate_output(vcf_file, cytoband_info, output_file)

