#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from pathlib import Path
from modules.params import initialize
from modules.readcount import launch_read_depth
from modules.normalize import launch_normalization
from modules.plot import plot_normalization, CnvPlot, plot_gene
from modules.cluster import launch_sample_clustering
from modules.ratio import calculate_coverage_ratios
from modules.segment import cbs, gaussian_hmm, custom_hmm_seg
from modules.call import call_cnvs, call_cnvs_2, export_all_calls
from modules.hmm import calculate_positional_mean_variance, CustomHMM

import numpy as np
import pandas as pd

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(main_dir,"/modules"))

def main(args):

    # Create output directory
    if not os.path.isdir(args.output_dir):
        os.mkdir(os.path.normpath(args.output_dir))

    output_name = os.path.basename(os.path.normpath(args.output_dir))
    log_file_name = ("{}{}").format(output_name,".grapes2.log")
    log_file = str(Path(args.output_dir) / log_file_name)

    # logging formatting
    logging.basicConfig(filename=log_file, filemode='w',
        format='%(asctime)s\t%(message)s')
    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # I/O Initialization
    sample_list, analysis_dict, ngs_utils_dict, ann_dict = initialize(args)

    # Depth, gc, mappability extraction formatting
    sample_list, analysis_dict = launch_read_depth(sample_list, analysis_dict,
        ngs_utils_dict, ann_dict)

    # Normalize by gc and then by mappability
    sample_list, analysis_dict = launch_normalization(sample_list, analysis_dict, ann_dict)

    # Plot nromaliation profiles
    # sample_list = plot_normalization(sample_list, analysis_dict)

    sample_list, analysis_dict = launch_sample_clustering(sample_list, analysis_dict)

    # calculate ratios
    sample_list, analysis_dict = calculate_coverage_ratios(sample_list, analysis_dict)

    genes = analysis_dict['list_genes_to_plot']
    for gene in genes:
        for sample in sample_list:
            plot_gene(sample.name, sample_list, gene, analysis_dict)
    sys.exit()

    # obs_dict = calculate_positional_mean_variance(sample_list, analysis_dict)
    sample_list = custom_hmm_seg(sample_list, analysis_dict)

    sample_list = call_cnvs_2(sample_list, -0.621, 0.433)
    export_all_calls(sample_list, analysis_dict)
    # sample_list = call_cnvs(sample_list, -0.621, 0.433, 2.58)

    for sample in sample_list:
        # Instantiating a CnvPlot object
        cnp = CnvPlot(
            cnr_file  = sample.ratio_file,
            cns_file  = ".",
            calls     = ".",
            sample    = sample.name,
            output_dir= sample.sample_folder,
            dup_cutoff = 0.433,
            del_cutoff = -0.621
        )
        # Plotting genomewide CNV profile
        sample_plot = cnp.plot_genomewide(genomewide=True, by_chr=False)

def parse_arguments():

    parser = argparse.ArgumentParser(description='GRAPES2: Detection of CNVs on gene panels')
    parser.add_argument('--bam_dir', type=str, dest='bam_dir', required=True,
        help='Input directory with bam files to be analyzed')
    parser.add_argument('--output_dir', dest='output_dir', required=True,
        help='Output directory')
    parser.add_argument('--bed', dest='bed', required=True,
        help='BED target regions')
    parser.add_argument("-t", "--threads", type=int, default=4,
        help="Number of CPU threads", dest='threads')
    parser.add_argument("--database", type=str, dest='database',
        help="Database directory harbouring sqlite files")
    parser.add_argument("-f", "--fasta", required=True, type=str,
        help="Genome reference in FASTA format", dest='reference')
    parser.add_argument("--upper_del_cutoff", type=float, default=-0.5,
        help=".", dest='upper_del_cutoff')
    parser.add_argument("--lower_del_cutoff", type=float, default=-1.5,
        help=".", dest='lower_del_cutoff')
    parser.add_argument('--plot_normalization', default=False,
        help="Plot depth normalization")
    parser.add_argument('--plot_gene', default=None,
        help="Plot gene log2 ratios by exon", dest='plot_gene')
    # parser_all.add_argument('--breakpoint', default=True,
    #     help="Perform offtarget breakpoint")

    args = parser.parse_args()
    return args

#########################################################
if __name__ == '__main__':
    args = parse_arguments()
    main(args)
