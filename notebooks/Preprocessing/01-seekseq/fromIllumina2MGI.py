### Acceleration by parallel reading
#!/usr/bin/env python
# Python 3.7 
import gzip
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(prog = "fromIllumina2MGI")
    parser.add_argument("-i", "--input", required = True, help = "Input fastq <fastq1>,<fastq2>")
    parser.add_argument("-o", "--output", required = True, help = "Output File <out1>,<out2>")
    args = parser.parse_args()
    return args

def rename_reads(input_file, output_file, read_type):
    """
    Correctly rename the reads in a FASTQ file according to the provided format.

    read_type (str): '/1' for R1 or '/2' for R2 reads.
    """
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        for line in infile:
            if line.startswith('@'):
                # Extract the relevant parts from the read name and format correctly
                line = line.rstrip('\n')
                parts = line.split(':')
                new_name = f"@{parts[0][1:]}{parts[1]}{parts[2]}{parts[3]}{parts[4]}{parts[5]}{parts[6].replace(' ', '')}{read_type}\n"
                outfile.write(new_name)
            else:
                # Write other lines as is
                outfile.write(line)

if __name__ == "__main__":

    args = parseArgs()

    #input /home/zhanghr/seq/Direct-Circ/ex_A180-sc-circ1-circ2-2-direct
    #output ../data/A180-direct-circ
    input_1, input_2 = args.input.split(',')
    output_1, output_2 = args.output.split(',')
    rename_reads(input_1, output_1, '/1')
    rename_reads(input_2, output_2, '/2')