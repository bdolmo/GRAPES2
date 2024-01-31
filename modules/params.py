import os
import sys
from modules.utils import get_bam_files, check_executable
from modules.sample import Sample
from pathlib import Path
import logging
from collections import defaultdict
from datetime import datetime

main_dir = Path(__file__).resolve().parents[1]
bin_dir = os.path.join(main_dir, "bin")
ann_dir = os.path.join(main_dir, "annotations")


def load_ngs_utils_config():
    """ """
    ngs_utils_dict = {
        # "megadepth": os.path.join(bin_dir, "megadepth"),
        # "mosdepth": os.path.join(bin_dir, "mosdepth"),
        "grapes_sv": os.path.join(bin_dir, "grapes_sv", "GRAPES"),
        "targetdepth": os.path.join(bin_dir, "TargetDepth", "target_depth.py"),
        # "targetdepth": os.path.join(bin_dir, "TargetDepth", "targetDepth.pl"),
    }
    for soft in ngs_utils_dict:
        check_executable(ngs_utils_dict[soft], dump_messages=True)

    return ngs_utils_dict


def load_annotation_config(genome_version: str):
    """ """
    ann_dict = {
        "mappability": os.path.join(
            ann_dir, "mappability", "wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz"
        ),
        "blacklist": os.path.join(
            ann_dir, "blacklist", "consensusBlacklist.hg19.bed"
        ),
        "chromosomes": os.path.join(
            ann_dir, "chromosomes", "hg19.chromosomes.txt"
        )
    }
    if genome_version == "hg38":
        ann_dict = {
            "mappability": os.path.join(
                ann_dir, "mappability", "wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz"
            ),
            "blacklist": os.path.join(
                ann_dir, "blacklist", "consensusBlacklist.hg38.bed"
            ),
            "chromosomes": os.path.join(
                ann_dir, "chromosomes", "hg38.chromosomes.txt"
            )
        }

    for item in ann_dict:
        if not os.path.isfile(ann_dict[item]):
            raise FileNotFoundError(ann_dict[item])
        else:
            msg = f" INFO: Found annotation {item}"
            logging.info(msg)
    return ann_dict


def initialize(args):
    """ """
    ngs_utils_dict = load_ngs_utils_config()
    annotation_dict = load_annotation_config(args.genome_version)
    analysis_dict = vars(args)

    analysis_dict["output_name"] = os.path.basename(analysis_dict["output_dir"])
    analysis_dict["list_genes_to_plot"] = []
    if args.plot_gene is not None:
        analysis_dict["list_genes_to_plot"] = args.plot_gene.replace(" ", "").split(",")
    analysis_dict["force"] = args.force

    # load bam files from input directory
    bam_list = []
    if os.path.isdir(args.bam_dir):
        bam_list = get_bam_files(args.bam_dir)
    else:
        with open(args.bam_dir) as f:
            for line in f:
                line = line.rstrip("\n")
                bam_list.append(line)
        f.close()

    now = datetime.now()
    date_time = now.strftime("%Y%m%d")

    # Iterate through bam list and create sample objects
    sample_list = list()
    for bam_file in sorted(bam_list):
        # Sample name from bam file prefix. todo: get name from bam SN
        sample_name = os.path.basename(bam_file).replace(".bam", "")

        analysis_json = defaultdict(dict)
        analysis_json["analysis_date"] = date_time
        analysis_json["sample_name"] = sample_name

        # Now create a sample object
        sample = Sample(sample_name)
        sample_folder = Path(args.output_dir) / sample_name
        sample_folder.mkdir(parents=True, exist_ok=True)
        sample.add("sample_folder", str(sample_folder))
        sample.add("bam", bam_file)
        sample.add("analysis_json", analysis_json)
        sample_list.append(sample)

    return sample_list, analysis_dict, ngs_utils_dict, annotation_dict
