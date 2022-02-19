import os
import sys
from modules.utils import get_bam_files, check_executable
from modules.sample import Sample
from pathlib import Path
import logging

main_dir = Path(__file__).resolve().parents[1]
bin_dir  = os.path.join(main_dir, "bin")
ann_dir  = os.path.join(main_dir, "annotations")

def load_ngs_utils_config():
    '''
    '''
    ngs_utils_dict = {
        'megadepth' : os.path.join(bin_dir, "megadepth"),
        'mosdepth'  : os.path.join(bin_dir, "mosdepth"),
        'grapes_sv' : os.path.join(bin_dir, "grapes_sv", "GRAPES"),
        'targetdepth': os.path.join(bin_dir, "TargetDepth", "targetDepth.pl")
    }
    for soft in ngs_utils_dict:
        check_executable(ngs_utils_dict[soft], dump_messages=True)

    return ngs_utils_dict

def load_annotation_config():
    '''
    '''
    ann_dict = {
        'mappability': os.path.join(ann_dir, "mappability",
            "wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz")
    }
    for item in ann_dict:
        if not os.path.isfile(ann_dict[item]):
            raise FileNotFoundError(ann_dict[item])
        else:
            msg = (" INFO: Found annotation {}").format(item)
            logging.info(msg)
    return ann_dict

def initialize(args):
  '''
  '''
  ngs_utils_dict  = load_ngs_utils_config()
  annotation_dict = load_annotation_config()
  analysis_dict   = vars(args)

  analysis_dict['output_name'] = os.path.basename(analysis_dict['output_dir'])

  analysis_dict['list_genes_to_plot'] = args.plot_gene.replace(" ", "").split(",")


  # load bam files from input directory
  bam_list = get_bam_files(args.bam_dir)

  # Iterate through bam list and create sample objects
  sample_list = list()
  for bam_file in sorted(bam_list):
    # Sample name from bam file prefix. todo: get name from bam SN
    sample_name = os.path.basename(bam_file).replace(".bam", "")

    # Now create a sample object
    sample = Sample(sample_name)
    sample_folder = Path(args.output_dir) / sample_name
    sample_folder.mkdir(parents=True, exist_ok=True)
    sample.add("sample_folder", str(sample_folder))
    sample.add("bam", bam_file)
    sample_list.append(sample)

  return sample_list, analysis_dict, ngs_utils_dict, annotation_dict
