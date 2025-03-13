import Levenshtein as lv
import sys
import time
import gzip
import argparse
from datetime import datetime

class Logger:
    def __init__(self, log_file_name, stream=sys.stdout):
        self.log_file_name = log_file_name
        self.stream = stream

    def write(self, message):
        with open(self.log_file_name, 'a') as log_file:
            log_file.write(message)
        if self.stream:
            self.stream.write(message)

    def flush(self):
        if self.stream:
            self.stream.flush()

########
# change in 2023-12-30
# More elegant
# Add log file
# Add event_counts_file  
def extract_paired_sequences(r1_input_file, fq1_output_file, r2_input_file, fq2_output_file, const="TTGTCTTCCTAAGAC", scaffold="GTTTTAGAGCTAGAA", max_dist=1):
    # Logger file generation
    log_file_name = f"seek_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    logger = Logger(log_file_name)

    len_const = len(const)
    len_sc = len(scaffold)

    event_counts = []
    current_count = 0

    with gzip.open(r1_input_file, 'rt') as r1_file, gzip.open(fq1_output_file, 'wt') as fq1_output, \
         gzip.open(r2_input_file, 'rt') as r2_file, gzip.open(fq2_output_file, 'wt') as fq2_output:

        read_count = 0
        matched_pairs = 0

        while True:
            read_count += 1

            # Read in R1 R2
            r1_identifier = r1_file.readline().strip()
            r1_sequence = r1_file.readline().strip()
            r1_plus_line = r1_file.readline().strip()
            r1_quality = r1_file.readline().strip()

            r2_identifier = r2_file.readline().strip()
            r2_sequence = r2_file.readline().strip()
            r2_plus_line = r2_file.readline().strip()
            r2_quality = r2_file.readline().strip()

            if not (r1_quality and r2_quality):  # End of file
                break
            
            # Record UMI
            if read_count % 10000 == 0:
                event_counts.append(current_count)
                logger.write(f"Processed {read_count} pairs of reads... Currently hits {matched_pairs}...\n")
            
            r1_match = False
            r2_match = False

            # Check R1
            for i in range(len(r1_sequence) - len_const + 1):
                if lv.distance(const, r1_sequence[i:i+len_const]) <= max_dist:
                    r1_start_index = max(0, i-25)
                    r1_end_index = min(len(r1_sequence), i + len_const + 10)
                    event_sequence = r1_sequence[i+len_const:r1_end_index]
                    r1_extracted_seq = r1_sequence[r1_start_index:i] + r1_sequence[i+len_const:r1_end_index]
                    r1_extracted_quality = r1_quality[r1_start_index:i] + r1_quality[i+len_const:r1_end_index]
                    if len(r1_extracted_seq) == 35:
                        r1_match = True
                    break

            # Check R2
            for i in range(len(r2_sequence) - len_sc + 1, 0, -1): # Search from end would be faster
                if lv.distance(scaffold, r2_sequence[i:i+len_sc]) <= max_dist and i-21 >= 0:
                    r2_extracted_seq = r2_sequence[i-21:i]
                    r2_extracted_quality = r2_quality[i-21:i]
                    r2_match = True
                    break
            
            if r1_match and r2_match:
                # # check length
                # # if len(r1_extracted_seq) == 35 and len(r2_extracted_seq) == 21: # redundant
                # is_new_event = True
                # for existing_event in unique_events:
                #     if lv.distance(event_sequence, existing_event) <= 1:
                #         is_new_event = False
                #         break
                # if is_new_event:
                #     unique_events.add(event_sequence)
                #     current_count += 1

                matched_pairs += 1
                fq1_output.write(f"{r1_identifier}\n{r1_extracted_seq}\n{r1_plus_line}\n{r1_extracted_quality}\n")
                fq2_output.write(f"{r2_identifier}\n{r2_extracted_seq}\n{r2_plus_line}\n{r2_extracted_quality}\n")

        #event_counts_file_name = f"event_counts_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        #with open(event_counts_file_name, 'w') as file:
        #    for count in event_counts:
        #       file.write(f"{count}\n")      

        logger.write(f"Finished processing. Total pairs of reads processed: {read_count}. Total matched pairs: {matched_pairs}\n")
        # logger.write(f"Total unique events: {current_count}\n")

# 使用示例
# extract_paired_sequences('r1_input.gz', 'fq1_output.txt', 'r2_input.gz', 'fq2_output.txt')
########
# change in 2023-12-07
# R1 is CID and MID region
# R2 is close to guide
def parseArgs():
    parser = argparse.ArgumentParser(prog = "SeekSeq")
    parser.add_argument("-i", "--input", type = str, required = True, help = "Input fastq File: <in1>,<in2>")
    parser.add_argument("-o", "--output", type = str, required = True, help = "Output fastq File: <out1>,<out2>")
    parser.add_argument("-c", "--constant", type = str, help = "Constant Sequence in read1", default = "TTGTCTTCCTAAGAC")
    parser.add_argument("-s", "--scaffold", type = str, help = "Scaffold Sequence in read2", default = "GTTTTAGA")
    parser.add_argument("-e", "--maxerror", type = int, help = "Max wrong bases to match constant and scaffold", default = 1)
    args = parser.parse_args()
    return(args)

if __name__ in "__main__":

    args = parseArgs()

    input_1, input_2 = args.input.split(',')
    output_1, output_2 = args.output.split(',')

    # example
    # extract_paired_sequences('./A180-sc-circ1-circ2_R2.fq.gz', './ex_A180-sc-circ1-circ2-2-direct_R1.fq', './A180-sc-circ1-circ2_R1.fq.gz', './ex_A180-sc-circ1-circ2-2-direct_R2.fq', const = "TTGTCTTCCTAAGAC", scaffold = "GTTTTAGA",max_dist =1)
    start_time = time.time()
    extract_paired_sequences(input_1, output_1, input_2, output_2, const = "TTGTCTTCCTAAGAC", scaffold = "GTTTTAGA",max_dist =1)
    end_time = time.time()
    print("Total execution time: {:.2f} seconds".format(end_time - start_time))