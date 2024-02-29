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

    if args.use_baseline_db:
        if not args.baseline_db:
            msg = "Missing baseline_db file (--baseline_db param must be filled)"
            print(msg)
            sys.exit()

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

    splitted_bed = args.bed.replace(".bed", ".splitted.bed")
    split_large_exons(args.bed, splitted_bed)
    analysis_dict["bed"] = splitted_bed


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


def split_large_exons(input_bed, output_bed):
    with open(input_bed, "r") as f:
        with open(output_bed, "w") as o:
            for line in f:
                # Split the line by tab character
                tmp = line.strip().split("\t")
                start = int(tmp[1])
                end = int(tmp[2])
                exon = tmp[3]
                size = end - start
                
                if size > 250:
                    # Calculate the number of segments required
                    num_segments = size // 250
                    
                    # Calculate the size of each segment
                    segment_size = size // num_segments
                    
                    # Write the segments to the output bed file
                    for i in range(num_segments):
                        segment_start = start + i * segment_size
                        segment_end = segment_start + segment_size
                        o.write("\t".join([tmp[0], str(segment_start), str(segment_end), exon]))
                        o.write("\n")
                        
                    # If there's a remaining part, write it as a separate segment
                    remaining_size = size % num_segments
                    if remaining_size > 0:
                        segment_start = start + num_segments * segment_size
                        segment_end = end
                        o.write("\t".join([tmp[0], str(segment_start), str(segment_end), exon]))
                        o.write("\n")
                else:
                    # If the exon size is not greater than 250 bp, just write it to the output file
                    o.write(line)