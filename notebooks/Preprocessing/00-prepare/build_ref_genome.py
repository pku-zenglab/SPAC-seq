#!/usr/bin/env python
# Python 3.7 
# without using U6 + scaffold
import argparse
import pandas as pd

# Global variables
fasta_lines = []
gtf_lines = []

# Taken from crop-seq
fasta_header_template = ">chr{chrom}\t{chrom}"

gtf_template = """chr{chrom}\tHAVANA\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_version "5"; gene_type "lncRNA"; gene_name "{id}_gene"; level 1; mgi_id "MGI:{id}"; havana_gene "{id}_gene";
chr{chrom}\tHAVANA\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_version "5";  transcript_id "{id}_transcript"; transcript_version "1"; gene_type "lncRNA"; gene_name "{id}_gene"; transcript_type "lncRNA"; transcript_name "{id}_transcript"; level 1; transcript_support_level "1"; mgi_id "MGI:{id}"; havana_gene "{id}_gene"; havana_transcript "{id}_transcript";
chr{chrom}\tHAVANA\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_version "5";  transcript_id "{id}_transcript"; transcript_version "1"; gene_type "lncRNA"; gene_name "{id}_gene"; transcript_type "lncRNA"; transcript_name "{id}_transcript"; exon_number 1; exon_id "{id}_exon"; exon_version "1"; level 1; transcript_support_level "1"; mgi_id "MGI:{id}"; havana_gene "{id}_gene"; havana_transcript "{id}_transcript";
"""

# Build sequence
def buildRef(row):
    # Get row with these columns from the CSV file of gRNAs
    sgRNA = row["sgRNA"]
    ref_seq = row["sequence"]
    
    # WRITE FASTA
    header = fasta_header_template.format(chrom=sgRNA)
    fasta_lines.append(header)
    fasta_lines.append(ref_seq)

    # WRITE GTF
    gtf_lines.append(gtf_template.format(chrom=sgRNA, id=sgRNA, length=len(ref_seq)))
    
# Parsers
def loadFiles(args):
    guide_df = pd.read_csv(args.input)
    return(guide_df)
# Argument parser
def parseArgs():
    parser = argparse.ArgumentParser(prog = "PrepRef")
    parser.add_argument("-i", "--input", type = str, required = True, help = "Input gRNA file")
    parser.add_argument("-o", "--output", type = str, required = True, help = "Output files name")
    args = parser.parse_args()
    return(args)

if __name__ ==  "__main__":
    # Parse arguments
    args = parseArgs()

    # Parse spreadsheets
    guide_df = loadFiles(args)

    guide_df.apply(buildRef, axis = 1)

    # Write to file
    output_fasta = args.output + ".fa"
    output_gtf = args.output + ".gtf"
    
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join(fasta_lines))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_lines)