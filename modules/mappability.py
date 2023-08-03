#!/usr/bin/env python3

import os
import sys
import re
import logging
import gzip
from datetime import datetime
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np
import subprocess
from natsort import natsorted, index_natsorted, order_by_index
import pybedtools
from collections import defaultdict


def annotate_mappability_bed(input_bed, mappability_bed):
    """ """
    output_bed = input_bed.replace(".bed", ".map.bed")
    if not os.path.isfile(output_bed):

        msg = " INFO: Extracting off-target mappability"
        logging.info(msg)

        a = pybedtools.BedTool(mappability_bed)
        b = pybedtools.BedTool(input_bed)
        c = a.intersect(b, wo=True, stream=True)
        map_dict = defaultdict(dict)

        for line in iter(c):
            line = str(line)
            line = line.rstrip()
            tmp = line.split("\t")
            mappability = float(tmp[3])
            size = int(tmp[6]) - int(tmp[5])
            bases = int(tmp[-1])
            marginal = ((mappability * bases) / size) * 100
            coordinate = "\t".join(tmp[4:9])
            if not coordinate in map_dict:
                map_dict[coordinate] = marginal
            else:
                map_dict[coordinate] += marginal

        o = open(output_bed, "w")
        for region in map_dict:
            o.write(region + "\t" + str(map_dict[region]) + "\n")
        o.close()
    return output_bed


def annotate_mappability(analysis_dict, ann_dict):
    """ """

    msg = " INFO: Extracting mappability"
    logging.info(msg)

    gc_bed = analysis_dict["gc_bed"]

    map_bed_name = os.path.basename(gc_bed).replace(".gc.bed", ".map.bed")
    map_bed = str(Path(analysis_dict["output_dir"]) / map_bed_name)

    if not os.path.isfile(map_bed):
        a = pybedtools.BedTool(ann_dict["mappability"])
        b = pybedtools.BedTool(gc_bed)
        c = a.intersect(b, wo=True, stream=True)
        map_dict = defaultdict(dict)

        for line in iter(c):
            line = str(line)
            line = line.rstrip()
            tmp = line.split("\t")
            mappability = float(tmp[3])
            size = int(tmp[6]) - int(tmp[5])
            bases = int(tmp[-1])
            marginal = ((mappability * bases) / size) * 100
            coordinate = "\t".join(tmp[4:9])
            if not coordinate in map_dict:
                map_dict[coordinate] = marginal
            else:
                map_dict[coordinate] += marginal

        o = open(map_bed, "w")
        for region in map_dict:
            o.write(region + "\t" + str(map_dict[region]) + "\n")
        o.close()
    else:
        msg = " INFO: Skipping mappability extraction"
        logging.info(msg)

    analysis_dict["map_bed"] = map_bed
    analysis_dict["ready_bed"] = map_bed

    return analysis_dict
