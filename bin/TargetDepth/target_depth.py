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
from collections import defaultdict

def which(program):
    return shutil.which(program)

def run_cmd(command):
    subprocess.run(command, shell=True, 
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def get_total_reads(bam_file):
    total_reads = 0
    for stat in pysam.idxstats(bam_file).splitlines():
        fields = stat.split('\t')
        total_reads += int(fields[2])
    return total_reads


def process_bam(bam_args):
    bam_file, args, depth_option = bam_args
    """ """
    # Determine the directory in which the Python script is located
    dirname = os.path.dirname(os.path.abspath(__file__))

    bam_name = os.path.basename(bam_file)

    sample_counts_file = f"{args.outdir}/{bam_name}_counts.bed"
    sample_coverage_file = f"{args.outdir}/{bam_name}_coverage.bed"
    # print(sample_counts_file, sample_coverage_file)

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


def merge_files(file_list, output_file):
    """ """
    # Extract sample names from file names
    sample_names = [os.path.basename(f).replace("_counts.bed", "")\
        .replace("_coverage.bed", "").replace(".bam", "") for f in file_list]

    # Create header with sample names
    header = "chr\tstart\tend\texon\tgc\tmap\t" + "\t".join(sample_names)
    sample_depth_dict = {}
    sample_depth = []

    for input_file in file_list:
        sample_name = os.path.basename(input_file).replace("_counts.bed", "").replace(".bam", "")

        if not sample_name in sample_depth_dict:
            sample_depth_dict[sample_name] = []

        with open (input_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("chr\t"):
                    continue
                tmp = line.split("\t")
                depth_data = tmp[-1]
                sample_depth_dict[sample_name].append(depth_data)
                # for data in sample_depth_dict:
                #     sample_depth_dict[sample_name].append(depth_data)

        f.close()
    
    baseline_file = file_list[0]
    o = open(output_file, "w")
    o.write(header+"\n")
    idx = 0
    with open(baseline_file) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            list_depth = []
            for sample in sample_depth_dict:
                sample_depth = sample_depth_dict[sample][idx]
                list_depth.append(sample_depth)
            o.write('\t'.join(tmp[0:6])+ "\t" + '\t'.join(list_depth)+"\n")
            idx+=1
    o.close()

    # print(file_list)


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
        print(f" ERROR: missing {summary_file} file")
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
                sample = data[0]
                sample_bam = os.path.join(args.input, sample)

                total_reads = get_total_reads(sample_bam)

                reads_on_target = int(data[1])
                reads_X = int(data[2])
                mean_coverage = float(data[3])
                mean_counts = float(data[4])
                mean_isize = float(data[5]) 
                sd_isize = float(data[6])
                mean_coverage_X = float(data[7])
                mean_counts_X = float(data[8]) 

                roi = (reads_on_target / total_reads) * 100 if total_reads > 0 else 0

                sf.write(f"{sample}\t{total_reads}\t{reads_on_target}\t{reads_X}\t{roi}\t{mean_coverage}\t{mean_counts}\t{mean_isize}\t{sd_isize}\t{mean_coverage_X}\t{mean_counts_X}\n")
        rf.close()
        sf.close()
        os.remove(summary_file)
        os.rename(new_summary_file, summary_file)

    # Merging temporary files
    count_files = glob.glob(os.path.join(args.outdir, "*_counts.bed"))
    coverage_files = glob.glob(os.path.join(args.outdir, "*_coverage.bed"))

    master_counts = os.path.join(args.outdir, f"{args.output_name}.read.counts.bed")
    master_coverage = os.path.join(args.outdir, f"{args.output_name}.per.base.coverage.bed")

    # Assuming that the first file in each list contains the header
    # Modify as needed based on your file format
    if args.report_counts:
        merge_files(count_files, master_counts)
    if args.report_coverage:
        merge_files(coverage_files, master_coverage)

    # Clean up temporary files
    # for f in count_files + coverage_files:
    #     os.remove(f)

    # temp_files = glob.glob(os.path.join(args.outdir, "*.tmp"))
    # for temp_file in temp_files:
    #     os.remove(temp_file)

    # Final message or summary
    print("Processing complete. Results are available in:", args.outdir)


if __name__ == "__main__":
    main()
