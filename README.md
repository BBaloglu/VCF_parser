# VCF_parser
Generates chromosome arm information using VCF and cytoband files

Dependencies:
Requires installation of pyvcf, A Variant Call Format Parser for Python

Install as follows:
pip3 install pycvf

See requirements.txt for the other packages used in the code

Run the script as follows:
python3 main.py [-i INFILE [INFILE ...]] [-c CYTOBAND] [-o OUTFILE]

Example command:
python3 main.py -i CN_segments.vcf -c ../cytoBand.txt -o cytoband_out.tsv
