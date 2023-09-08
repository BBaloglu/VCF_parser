#!/usr/bin/python3

import vcf

# Define the path to cytoband.txt (modify as needed)
cytoband_file = '../cytoband.txt'

def read_cytoband_info(cytoband_file):
    """
    Read cytoband information from the cytoband.txt file and store it in a dictionary.

    Args:
        cytoband_file (str): Path to the cytoband.txt file.

    Returns:
        dict: A dictionary containing cytoband information for each chromosome.
    """
    cytoband_info = {}
    with open(cytoband_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            band = parts[3]
            cytoband_info.setdefault(chromosome, []).append((start, end, band))
    return cytoband_info

def find_cytoband_range(chromosome, start, end, cytoband_info):
    """
    Find the cytoband range for a given chromosome, start, and end positions.

    Args:
        chromosome (str): Chromosome name (e.g., "chr1").
        start (int): Start position in the chromosome.
        end (int): End position in the chromosome.
        cytoband_info (dict): Dictionary containing cytoband information.

    Returns:
        str: Cytoband range (e.g., "p11.32q23").
    """
    if chromosome in cytoband_info:
        cytoband_ranges = []
        for s, e, band in cytoband_info[chromosome]:
            if s <= start <= e or s <= end <= e:
                cytoband_ranges.append(band)
        if cytoband_ranges:
            return ''.join(sorted(cytoband_ranges))
    return "Unknown"

def parse_cnv_and_generate_output(vcf_file, output_file, cytoband_file):
    """
    Compare chromosome information from a VCF file to cytoband information from a text file
    and generate an output file with chromosome arm band information and ploidy.

    Args:
        vcf_file (str): Path to the VCF file.
        output_file (str): Path to the output file.
        cytoband_file (str): Path to the cytoband.txt file.

    Returns:
        None
    """
    try:
        cytoband_info = read_cytoband_info(cytoband_file)

        # Open VCF file using PyVCF
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))

        # Open the output file for writing
        with open(output_file, 'w') as output:
            # Write the header for the output file
            output.write("ChrNbr\tLength\tPloidy\tCytoband\n")

            # Process each record in the VCF file
            for record in vcf_reader:
                chromosome = record.CHROM
                start_position = record.POS
                print('start_pos', start_position)
                end_position = record.INFO['END']
                print('end position', end_position)
                length = end_position - start_position
                print('length', length)
                ploidy = record.genotype(record.samples[0].sample)['CN']  # Extract ploidy from the first sample
                cytoband_range = find_cytoband_range(chromosome, start_position, end_position, cytoband_info)
                print('cytoband_range', cytoband_range)

                # Write the information to the output file
                output.write(f"{chromosome}\t{length}\t{ploidy:.2f}\t{chromosome}{cytoband_range}\n")

        print("Output file generated successfully.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Prompt the user to enter the path to the VCF file
    vcf_file = input("Enter the path to the VCF file: ")
    output_file = input("Enter the path for the output file: ")

    print("Parsing cytoband information...")
    parse_cnv_and_generate_output(vcf_file, output_file, cytoband_file)

