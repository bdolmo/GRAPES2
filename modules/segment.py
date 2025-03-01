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
import subprocess
from hmmlearn import hmm
from scipy.special import logsumexp


#from modules.hmm import calculate_positional_mean_variance, CustomHMM
from modules.hmmnb import calculate_positional_mean_variance, CustomHMM


def custom_hmm_seg(sample_list, analysis_dict):
    """ """
    obs_dict = calculate_positional_mean_variance(sample_list, analysis_dict)

    for sample in sample_list:

        if sample.analyzable == "False":
            continue

        chr_dict = load_observations_by_chr(sample.ratio_file)

        segment_file_name = f"{sample.name}.segment.bed"
        segment_file = str(Path(sample.sample_folder) / segment_file_name)
        sample.add("segment_file", segment_file)

        segment_file_extended_name = f"{sample.name}.extended.segment.bed"
        segment_file_extended = str(
            Path(sample.sample_folder) / segment_file_extended_name
        )
        sample.add("segment_extended_file", segment_file_extended)

        segment_file_map_name = f"{sample.name}.segment.map.bed"
        segment_file_map = str(
            Path(sample.sample_folder) / segment_file_map_name
        )
        sample.add("segment_file_map", segment_file_map)
        if os.path.isfile(segment_file):
            continue

        m = open(segment_file_map, "w")
        o = open(segment_file, "w")
        p = open(segment_file_extended, "w")

        msg = f" INFO: segmenting sample {sample.name}"
        logging.info(msg)

        for chr in chr_dict:

            model = CustomHMM(obs_dict, sample.name, chr)

            dispersions = model.fit_dispersion(max_iter=10, tol=1e-3)
            alpha = model.forward()
            print(sample.name, chr, "dispersion_values", dispersions)

            states, phred_scores = model.decode()
            posteriors = model.posterior_decoding()

            # Calculate MAP estimates
            map_states, map_probabilities = model.calculate_map()
            log_likelihoods = model.compute_log_likelihood()
            
            # Output results
            jdx = 0
            for i, (state_map, prob_map) in enumerate(zip(map_states, map_probabilities)):
                # state_map = state_map.rstrip("\n")
                item = chr_dict[chr][jdx]
                posterior = str(posteriors[jdx]).rstrip("\n")
                tmp = item["region"].split("\t")

                # chr15	48704765	48704940	NM_000138_64_65;FBN1	52.0	100.0	-0.813	[-4.040678, 1.361648, -2.714965, -16.251632, -18.420681]	4	0.6206
                # Convert log likelihoods to a probability vector:
                logL_vector = np.array(log_likelihoods[jdx])
                prob_vector = np.exp(logL_vector - logsumexp(logL_vector))
                most_state = np.argmax(prob_vector)
                most_prob = prob_vector[most_state]
                
                m.write(f"{item['region']}\t{str(log_likelihoods[jdx])}\t{most_state}\t{most_prob:.4f}"+"\n")
                jdx += 1
            idx = 0
            unmerged_list = []
            for item in chr_dict[chr]:
                print(item, states, idx)
                state = states[idx]
                phred = phred_scores[idx][int(state)]              
                tmp = item["region"].split("\t")
                data_dict = {
                    "chr": tmp[0],
                    "start": tmp[1],
                    "end": tmp[2],
                    "region": tmp[3],
                    "gc": tmp[4],
                    "map": posteriors[idx],
                    "log2_ratio": tmp[6],
                    "state": str(state),
                    "phred": phred
                }
                unmerged_list.append(data_dict)
                p.write(item["region"] + "\t" + str(state) + "\t" + str(map_probabilities[idx])+ "\t" + str(state) + "\t" + str(posteriors[idx]) + "\n")
                idx += 1
            merged_list = merge_segments(unmerged_list)
            for item in merged_list:
                out_list = []
                for val in item:
                    out_list.append(str(item[val]))
                o.write("\t".join(out_list) + "\n")  
        o.close()
        p.close()

    return sample_list

def merge_segments(unmerged_list):
    """ """
    merged_list = []
    merging_items = []

    flag = 0
    for region in unmerged_list:
        if flag == 0:
            first_dict = region
            flag = 1
            merging_items.append(first_dict)
            continue
        second_dict = region
        if first_dict["state"] == second_dict["state"]:
            if first_dict is not second_dict:
                merging_items.append(second_dict)
        else:
            min_start = 10e20
            max_end = 0
            region_list = []
            ratio_list = []
            phred_list = []
            for item in merging_items:
                if int(item["start"]) < min_start:
                    min_start = int(item["start"])
                if int(item["end"]) > max_end:
                    max_end = int(item["end"])
                ratio_list.append(float(item["log2_ratio"]))
                region_list.append(item["region"])
                phred_list.append(item["phred"])

            # Compute average posterior probability
            mean_phred = np.mean(phred_list)  
            mean_ratio = round(np.median(ratio_list), 3)
            new_segment = {
                "chr": first_dict["chr"],
                "start": min_start,
                "end": max_end,
                "region": ",".join(region_list),
                "n_regions": str(len(region_list)),
                "log2_ratio": mean_ratio,
                "state": first_dict["state"],
                "phred": mean_phred
            }
            merged_list.append(new_segment)
            first_dict = region
            merging_items = []
            merging_items.append(second_dict)

    if len(merging_items) > 0:

        min_start = 10e20
        max_end = 0
        region_list = []
        ratio_list = []
        phred_list = []

        for item in merging_items:
            if int(item["start"]) < min_start:
                min_start = int(item["start"])
            if int(item["end"]) > max_end:
                max_end = int(item["end"])

            ratio_list.append(float(item["log2_ratio"]))
            region_list.append(item["region"])
            phred_list.append(item["phred"])

        mean_phred = np.mean(phred_list)  
        mean_ratio = round(np.median(ratio_list), 3)
        new_segment = {
            "chr": first_dict["chr"],
            "start": min_start,
            "end": max_end,
            "region": ",".join(region_list),
            "n_regions": str(len(region_list)),
            "log2_ratio": mean_ratio,
            "state": first_dict["state"],
            "phred": mean_phred
        }
        merged_list.append(new_segment)

    return merged_list


def load_observations_by_chr(ratio_file):
    """ """
    # Load observations
    chr_obs_dict = defaultdict(dict)
    with open(ratio_file) as f:
        for line in f:
            if line.startswith("chr\tstart"):
                continue
            line = line.rstrip("\n")
            tmp = line.split("\t")
            coordinate = ("{}\t{}\t{}\t{}\t{}\t{}\t{}\n").format(
                tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[-1]
            )
            chr = tmp[0]
            log2_ratio = tmp[-1]
            if not chr in chr_obs_dict:
                chr_obs_dict[chr] = []
            data_dict = {"log2_ratio": log2_ratio, "region": coordinate.rstrip("\n")}
            chr_obs_dict[chr].append(data_dict)

        f.close()
    return chr_obs_dict


def cbs(sample_list, n_segments=2, alpha=0.05):
    """
    Segment with CBS
    """
    for sample in sample_list:

        to_segment = sample.ratio_file.replace(".ratios.bed", ".tosegment.bed")
        o = open(to_segment, "w")
        with open(sample.ratio_file) as f:
            for line in f:
                if line.startswith("chr\tstart"):
                    continue
                line = line.rstrip("\n")
                tmp = line.split("\t")
                outline = ("{}\t{}\t{}\t{}\t{}\n").format(
                    tmp[0], tmp[1], tmp[2], tmp[3], tmp[-1]
                )
                o.write(outline)
            # f.close()
        o.close()

        rscript = sample.ratio_file.replace(".ratios.bed", ".CBS.R")
        segment_file = sample.ratio_file.replace(".bed", ".segment.bed")
        sample.add("segment_file", segment_file)

        r = open(rscript, "w")
        r.write("library(DNAcopy)" + "\n")
        line = 'cn <- read.table("{}", header=F)'.format(to_segment)
        r.write(line + "\n")
        line = "CNA.object <-CNA( genomdat = cn[,5], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')"
        r.write(line + "\n")
        line = "CNA.smoothed <- smooth.CNA(CNA.object)"
        r.write(line + "\n")
        line = (
            "segs <- segment(CNA.object, verbose=0, min.width={}, alpha = {})".format(
                n_segments, alpha
            )
        )
        r.write(line + "\n")
        line = "segs2=segs$output"
        r.write(line + "\n")
        line = 'write.table(segs2[,2:6], file="{}",row.names=F, col.names=F, quote=F, sep="\t")'.format(
            segment_file
        )
        r.write(line + "\n")
        r.close()

        if not os.path.isfile(segment_file):

            msg = (" INFO: Segmenting sample {}").format(sample.name)
            logging.info(msg)

            cmd = ("Rscript {}").format(rscript)
            p1 = subprocess.run(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            output = p1.stdout.decode("UTF-8")
            error = p1.stderr.decode("UTF-8")

        # add header to segment file
        idx = 0
        tmp_segment = segment_file.replace(".bed", ".tmp.bed")
        o = open(tmp_segment, "w")
        with open(segment_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if idx == 0:
                    if not line.startswith("chr"):
                        o.write("chr\tstart\tend\tsegments\tratio\n")
                o.write(line + "\n")
                idx += 1
        f.close()
        o.close()

        os.remove(segment_file)
        os.remove(rscript)
        os.rename(tmp_segment, segment_file)

    return sample_list
