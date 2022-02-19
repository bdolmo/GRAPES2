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
import pybedtools


def annotate_gc(analysis_dict):
    '''
        Add gc content
    '''
    bed = analysis_dict['bed']

    gc_bed_name = os.path.basename(bed).replace(".bed", ".gc.bed")
    gc_bed = str(Path(analysis_dict['output_dir']) / gc_bed_name)

    msg = " INFO: Extracting gc content for roi file {}".format(bed)
    a = pybedtools.BedTool(bed)
    b = a.nucleotide_content(fi=analysis_dict['reference'], pattern="CG", C=True)

    o = open(gc_bed, 'w')
    for line in iter(b):
        line = str(line)
        line = line.rstrip()
        tmp  = line.split('\t')
        tmp[5] = str(round(100*float(tmp[5]), 2))
        indices = [0,1,2,3,5]
        new_line = "\t".join([tmp[i] for i in indices])
        o.write(new_line + '\n')
    o.close()

    analysis_dict['gc_bed'] = gc_bed
    return analysis_dict
