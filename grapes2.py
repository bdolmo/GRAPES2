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
from modules.breakpoint import call_structural_variants
from modules.merge_cnv_sv import merge_bed_files
from modules.offtarget import create_offtarget_bed, create_pseudowindows, extract_offtarget
# from modules.breakpoint import call_structural_variants, export_sv_calls
# from modules.blat import Blat

from modules.call import (
    call_cnvs,
    call_raw_cnvs,
    filter_single_exon_cnv,
    unify_raw_calls,
    export_cnv_calls,
)

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(main_dir, "/modules"))

def main(args):

    # Create output directory
    if not os.path.isdir(args.output_dir):
        os.mkdir(os.path.normpath(args.output_dir))

    output_name = os.path.basename(os.path.normpath(args.output_dir))
    log_file_name = f"{output_name}.grapes.log"
    log_file = str(Path(args.output_dir) / log_file_name)

    # logging formatting
    logging.basicConfig(
        filename=log_file, filemode="w", format="%(asctime)s\t%(message)s"
    )
    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # I/O Initialization
    sample_list, analysis_dict, ngs_utils_dict, ann_dict = initialize(args)

    if args.offtarget:
        analysis_dict = create_offtarget_bed(args.bed, args.output_dir, args.reference, 
            analysis_dict, ann_dict["mappability"], ann_dict["chromosomes"], ann_dict["blacklist"])

        extract_offtarget(sample_list, args.reference, ngs_utils_dict, analysis_dict)

    if args.breakpoint:
        for sample in sample_list:
            # if sample.name != "sample49.simulated":
            #     continue

            msg = f" INFO: Calling breakpoints on sample {sample.name}"
            logging.info(msg)

            call_structural_variants(sample.bam, args.bed, args.reference, args.output_dir, 
                sample, analysis_dict, ngs_utils_dict, ann_dict)
        # sys.exit()

    # Depth, gc, mappability extraction formatting
    sample_list, analysis_dict = launch_read_depth(
        sample_list, analysis_dict, ngs_utils_dict, ann_dict
    )

    # Normalize by gc and then by mappability
    sample_list, analysis_dict = launch_normalization(
        sample_list, analysis_dict, ann_dict
    )

    # Plot normalization profiles
    # sample_list = plot_normalization(sample_list, analysis_dict)

    sample_list, analysis_dict = launch_sample_clustering(sample_list, analysis_dict)

    # calculate ratios
    sample_list, analysis_dict = calculate_coverage_ratios(sample_list, analysis_dict)

    sample_list = custom_hmm_seg(sample_list, analysis_dict)

    sample_list = call_raw_cnvs(
        sample_list, analysis_dict, args.upper_del_cutoff, args.lower_dup_cutoff
    )
    
    filter_single_exon_cnv(sample_list, args.upper_del_cutoff, 
        args.lower_dup_cutoff, analysis_dict
    )

    sample_list = unify_raw_calls(sample_list)

    sample_list = call_cnvs(sample_list, -0.621, 0.433, 2.58)
    sample_list = export_cnv_calls(sample_list, analysis_dict)

    for sample in sample_list:
        cnp = CnvPlot(
            cnr_file=sample.ratio_file,
            cns_file=".",
            calls=".",
            sample=sample.name,
            output_dir=sample.sample_folder,
            dup_cutoff=args.lower_dup_cutoff,
            del_cutoff=args.upper_del_cutoff,
        )
        # Plotting genomewide CNV profile
        sample_plot = cnp.plot_genomewide(genomewide=True, by_chr=False)


    for sample in sample_list:

        sv_bed = os.path.join(args.output_dir, sample.name, 
            f"{sample.name}.GRAPES2.breakpoints.bed")
        cnv_bed = os.path.join(args.output_dir, sample.name, 
            f"{sample.name}.GRAPES2.cnv.bed")
        merged_bed = os.path.join(args.output_dir, f"{sample.name}.GRAPES2.bed")
        merge_bed_files(cnv_bed, sv_bed, merged_bed)



def parse_arguments():

    parser = argparse.ArgumentParser(
        description="GRAPES2: Detection of CNVs on gene panels"
    )
    parser.add_argument(
        "--bam_dir",
        type=str,
        dest="bam_dir",
        required=True,
        help="Input directory with bam files to be analyzed",
    )
    parser.add_argument(
        "--output_dir", 
        dest="output_dir", 
        required=True, 
        help="Output directory"
    )
    parser.add_argument(
        "--bed", 
        dest="bed", 
        required=True, 
        help="BED target regions")
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads",
        dest="threads",
    )
    parser.add_argument(
        "--database",
        type=str,
        dest="database",
        help="Database directory harbouring sqlite files",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        type=str,
        help="Genome reference in FASTA format",
        dest="reference",
    )
    parser.add_argument(
        "--breakpoint", 
        action="store_true", 
        help="Perform breakpoint (SV) analysis",
        dest="breakpoint"
    )
    parser.add_argument(
        "--offtarget", 
        action="store_true", 
        help="Perform offtarget CNV analysis",
        dest="offtarget"
    )
    parser.add_argument(
        "--force",
        action="store_true", 
        dest="force",
    )
    parser.add_argument(
        "--upper_del_cutoff",
        type=float,
        default=-0.621,
        help=".",
        dest="upper_del_cutoff",
    )
    parser.add_argument(
        "--lower_dup_cutoff",
        type=float,
        default=0.433,
        help=".",
        dest="lower_dup_cutoff",
    )
    parser.add_argument(
        "--plot_normalization", 
        default=False, 
        help="Plot depth normalization"
    )
    parser.add_argument(
        "--plot_gene",
        default=None,
        help="Plot gene log2 ratios by exon",
        dest="plot_gene",
    )


    args = parser.parse_args()
    return args


#########################################################
if __name__ == "__main__":
    args = parse_arguments()
    main(args)
