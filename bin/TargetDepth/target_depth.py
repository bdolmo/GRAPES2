#!/usr/bin/env python3

import argparse
import os
import subprocess
import shutil
import pandas as pd
import glob
import sys
import pysam
from multiprocessing import Pool
from pathlib import Path


def which(program):
    return shutil.which(program)

def run_cmd(command):
    subprocess.run(command, shell=True, 
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def get_total_reads(bam_file):
    total_reads = 0
    # print(bam_file)

    for stat in pysam.idxstats(bam_file).splitlines():
        # print(stat)
        fields = stat.split('\t')
        total_reads += int(fields[2])
        # total_reads += int(fields[2])
    # for stat in pysam.idxstats(bam_file):
    #     print("stat", stat)
    #     fields = stat.split('\t')
    #     total_reads += int(fields[2])
    return total_reads


def process_bam(bam_args):
    bam_file, args, depth_option = bam_args
    """ """
    # Determine the directory in which the Python script is located
    dirname = os.path.dirname(os.path.abspath(__file__))

    bam_name = os.path.basename(bam_file)

    sample_counts_file = f"{args.outdir}/{bam_name}_counts.bed"
    sample_coverage_file = f"{args.outdir}/{bam_name}_coverage.bed"
    print(sample_counts_file, sample_coverage_file)

    # Check for TargetDepth binary execution
    targetDepthExe = os.path.join(dirname, "TargetDepth")
    if not os.path.isfile(targetDepthExe) or not os.access(targetDepthExe, os.X_OK):
        print("ERROR: Cannot execute TargetDepth")
        sys.exit(1)

    # bam_file, args = bam_args
    print(f"INFO: Processing {bam_file}")

    # Assuming targetDepthExe is a callable script
    cmd = f"{targetDepthExe} -i {bam_file} {depth_option} -g {args.genome} -b {args.bed} -o {args.outdir}"
    print(f"INFO: CMD: {cmd}")
    # if args.verbose:
    #     print(f"INFO: CMD: {cmd}")
    if not os.path.isfile(sample_counts_file) and not os.path.isfile(sample_coverage_file):
        run_cmd(cmd)
    print(f"INFO: Finished {bam_file}")

def merge_files(file_list, output_file, header=None):
    with open(output_file, 'w') as outfile:
        if header:
            outfile.write(header + '\n')
        for fname in file_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Extract read depth metrics from a list of files")
    parser.add_argument('-i', '--input', required=True, help="Input BAM/s. Options: directory, comma-separated bams, or file with bam paths")
    parser.add_argument('-o', '--outdir', default='.', help="Output directory")
    parser.add_argument('-c', '--report_counts', action='store_true', help="Report exon counts")
    parser.add_argument('-d', '--report_coverage', action='store_true', help="Report per base coverage")
    parser.add_argument('-b', '--bed', required=True, help="BED regions file")
    parser.add_argument('-n', '--output_name', required=True, help="Output name")
    parser.add_argument('-g', '--genome', required=True, help="Genome file in FASTA format")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads")
    parser.add_argument('-v', '--version', action='store_true', help="Display version")
    parser.add_argument('--verbose', action='store_true', help="Verbose output")
    args = parser.parse_args()

    if args.version:
        print("targetDepth.py v1.0")
        exit()

    # Check BED consistency
    with open(args.bed, 'r') as bed_file:
        for line in bed_file:
            if len(line.strip().split('\t')) < 4:
                print("ERROR: bed file requires a fourth field containing region information (e.g gene name)")
                exit()

    # Determine depth option based on arguments
    depth_option = ''
    if args.report_counts and args.report_coverage:
        depth_option = '-c -d'
    elif args.report_counts:
        depth_option = '-c'
    elif args.report_coverage:
        depth_option = '-d'

    bams = []
    if os.path.isdir(args.input):
        print(f" INFO: Analyzing all bam files from {args.input} directory")
        bams = [str(bam_file) for bam_file in Path(args.input).glob("*.bam")]
    elif os.path.isfile(args.input):
        with open(args.input, 'r') as file:
            bams = [line.strip() for line in file if os.path.isfile(line.strip())]
    else:
        bams = args.input.split(',')

    if not bams:
        print(" ERROR: No valid input files found")
        exit()

    print(f"Number of threads: {args.threads}")

    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    bam_arguments = [(bam_file, args, depth_option) for bam_file in bams]

    print(f" INFO: Processing {len(bams)} BAM files with {args.threads} processes")
    with Pool(processes=args.threads) as pool:
        pool.map(process_bam, bam_arguments)

    # Generating summary metrics
    summary_file = os.path.join(args.outdir, "summary_metrics.log")
    new_summary_file = os.path.join(args.outdir, "summary_metrics.tmp.log")

    samtools_path = which("samtools")
    
    if not os.path.exists(summary_file) or os.path.getsize(summary_file) == 0:
        print(f"ERROR: missing {summary_file} file")
    else:
        with open(summary_file, 'r') as rf:
            lines = rf.readlines()

        with open(new_summary_file, 'w') as sf:

            if samtools_path:
                sf.write("SAMPLE\tTOTAL_READS\tREADS_ON_TARGET\tREADS_CHRX\t%ROI\tMEAN_COVERAGE\tMEAN_COUNTS\tMEAN_ISIZE\tSD_ISIZE\tMEAN_COVERAGE_X\tMEAN_COUNTS_X\n")
            else:
                sf.write("SAMPLE\tREADS_ON_TARGET\tREADS_CHRX\tMEAN_COVERAGE\tMEAN_COUNTS\tMEAN_ISIZE\tSD_ISIZE\tMEAN_COVERAGE_X\tMEAN_COUNTS_X\n")

            for line in lines:
                data = line.strip().split('\t')
                if line.startswith("SAMPLE\t"):
                    continue
                # Process each line as required
                sample = data[0]
                sample_bam = os.path.join(args.input, sample)

                total_reads = get_total_reads(sample_bam)


                reads_on_target = int(data[1])
                reads_X = int(data[2])
                mean_coverage = float(data[3])
                mean_counts = float(data[4])
                mean_isize = float(data[5])  # Additional field
                sd_isize = float(data[6])  # Additional field
                mean_coverage_X = float(data[7])  # Additional field
                mean_counts_X = float(data[8])  # Additional field

                # Add more fields as per your file format
                # For example, you might have additional metrics or calculations here
                roi = (reads_on_target / total_reads) * 100 if total_reads > 0 else 0

                # Write processed data to the new summary file
                sf.write(f"{sample}\t{total_reads}\t{reads_on_target}\t{reads_X}\t{roi}\t{mean_coverage}\t{mean_counts}\t{mean_isize}\t{sd_isize}\t{mean_coverage_X}\t{mean_counts_X}\n")
                # Include more fields as necessary
        rf.close()
        sf.close()
        os.remove(summary_file)
        os.rename(new_summary_file, summary_file)
        # sys.exit()

    # Merging temporary files
    count_files = glob.glob(os.path.join(args.outdir, "*_counts.bed"))
    coverage_files = glob.glob(os.path.join(args.outdir, "*_coverage.bed"))

    master_counts = os.path.join(args.outdir, f"{args.output_name}.read.counts.bed")
    master_coverage = os.path.join(args.outdir, f"{args.output_name}.per.base.coverage.bed")

    # Assuming that the first file in each list contains the header
    # Modify as needed based on your file format
    if args.report_counts:
        merge_files(count_files, master_counts, header="chr\tstart\tend\texon\tgc\tmap\tSAMPLES")
    if args.report_coverage:
        merge_files(coverage_files, master_coverage, header="chr\tstart\tend\texon\tgc\tmap\tSAMPLES")

    # Clean up temporary files
    for f in count_files + coverage_files:
        os.remove(f)

    temp_files = glob.glob(os.path.join(args.outdir, "*.tmp"))
    for temp_file in temp_files:
        os.remove(temp_file)

    # Final message or summary
    print("Processing complete. Results are available in:", args.outdir)

# def analyze_data(counts_file, coverage_file):
#     """
#     Example function for additional analysis on the counts and coverage data.
#     This could involve statistical analyses, data aggregation, etc.
#     """

#     # Load the data into pandas DataFrames
#     counts_df = pd.read_csv(counts_file, sep='\t')
#     coverage_df = pd.read_csv(coverage_file, sep='\t')

#     # Example analysis:
#     # Calculate mean coverage and counts for each region
#     counts_mean = counts_df.mean(axis=1)
#     coverage_mean = coverage_df.mean(axis=1)

#     # Add these as new columns to the DataFrame
#     counts_df['Mean_Count'] = counts_mean
#     coverage_df['Mean_Coverage'] = coverage_mean

#     # Save the updated dataframes back to new files
#     counts_df.to_csv(counts_file.replace('.bed', '.mean_counts.bed'), sep='\t', index=False)
#     coverage_df.to_csv(coverage_file.replace('.bed', '.mean_coverage.bed'), sep='\t', index=False)

#     # Further analysis can be added here as per your requirements.
#     # This could include more complex statistical tests, data visualizations,
    # or exporting data for use in other tools.


if __name__ == "__main__":
    main()
