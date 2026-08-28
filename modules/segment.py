#!/usr/bin/env python3

import os
import logging
from collections import defaultdict
from pathlib import Path
import numpy as np
import subprocess
from modules.utils import (
    atomic_output_path,
    stage_manifest_matches,
    update_stage_manifest,
    validate_delimited_file,
)

#from modules.hmm import calculate_positional_mean_variance, CustomHMM

# Testing with a negative binomial distribution
from modules.hmmnb import calculate_positional_mean_variance, CustomHMM


def custom_hmm_seg(sample_list, analysis_dict):
    """ """
    obs_dict = None

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
        segment_metadata = {
            "chromosomes": len(chr_dict),
            "confidence_schema": 2,
        }
        segment_outputs = [
            segment_file,
            segment_file_extended,
            segment_file_map,
        ]
        segment_outputs_valid = (
            not analysis_dict.get("force", False)
            and stage_manifest_matches(
                analysis_dict["output_dir"],
                f"segmentation:{sample.name}",
                segment_outputs,
                input_paths=[sample.ratio_file],
                metadata=segment_metadata,
            )
            and
            validate_delimited_file(
                segment_file,
                min_columns=8,
                min_data_rows=1,
                has_header=False,
            )
            and validate_delimited_file(
                segment_file_extended,
                min_columns=11,
                min_data_rows=1,
                has_header=False,
            )
            and validate_delimited_file(
                segment_file_map,
                min_columns=10,
                min_data_rows=1,
                has_header=False,
            )
        )
        if segment_outputs_valid:
            logging.info(
                f" INFO: Reusing validated segmentation outputs for {sample.name}"
            )
            continue

        msg = f" INFO: segmenting sample {sample.name}"
        logging.info(msg)
        if obs_dict is None:
            obs_dict = calculate_positional_mean_variance(sample_list, analysis_dict)

        with atomic_output_path(segment_file) as segment_tmp, atomic_output_path(
            segment_file_extended
        ) as extended_tmp, atomic_output_path(segment_file_map) as map_tmp:
            with open(segment_tmp, "w", encoding="utf-8") as segment_handle, open(
                extended_tmp, "w", encoding="utf-8"
            ) as extended_handle, open(map_tmp, "w", encoding="utf-8") as map_handle:
                for chromosome in chr_dict:
                    model = CustomHMM(obs_dict, sample.name, chromosome)

                    dispersions = model.fit_dispersion(max_iter=10, tol=1e-3)
                    model.forward()
                    states, _transition_scores = model.decode()
                    posteriors = model.posterior_decoding()

                    # Calculate MAP estimates
                    _, map_probabilities = model.calculate_map()
                    log_likelihoods = model.compute_log_likelihood()

                    state_posteriors = []
                    for jdx, item in enumerate(chr_dict[chromosome]):
                        state = int(states[jdx])
                        posterior_row = np.asarray(posteriors[jdx])
                        if posterior_row.ndim == 0:
                            state_posterior = float(posterior_row)
                        else:
                            state_posterior = float(posterior_row[state])
                        state_posterior = min(max(state_posterior, 0.0), 1.0)
                        state_posteriors.append(state_posterior)
                        map_handle.write(
                            f"{item['region']}\t{str(log_likelihoods[jdx])}\t"
                            f"{state}\t{state_posterior:.6f}\n"
                        )

                    unmerged_list = []
                    for idx, item in enumerate(chr_dict[chromosome]):
                        state = states[idx]
                        state_posterior = state_posteriors[idx]
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
                            "hmm_posterior": state_posterior,
                        }
                        unmerged_list.append(data_dict)
                        extended_handle.write(
                            item["region"]
                            + "\t"
                            + str(state)
                            + "\t"
                            + str(map_probabilities[idx])
                            + "\t"
                            + str(state)
                            + "\t"
                            + str(posteriors[idx])
                            + "\n"
                        )

                    merged_list = merge_segments(unmerged_list)
                    for item in merged_list:
                        segment_handle.write(
                            "\t".join(str(item[value]) for value in item) + "\n"
                        )

                    state_values, state_counts = np.unique(states, return_counts=True)
                    state_summary = ", ".join(
                        f"{int(state)}:{int(count)}"
                        for state, count in zip(state_values, state_counts)
                    )
                    logging.info(
                        " INFO: Segmented %s %s: targets=%d segments=%d "
                        "states={%s} dispersions=%s",
                        sample.name,
                        chromosome,
                        len(chr_dict[chromosome]),
                        len(merged_list),
                        state_summary,
                        np.asarray(dispersions).round(6).tolist(),
                    )

            if not (
                validate_delimited_file(
                    segment_tmp,
                    min_columns=8,
                    min_data_rows=1,
                    has_header=False,
                )
                and validate_delimited_file(
                    extended_tmp,
                    min_columns=11,
                    min_data_rows=1,
                    has_header=False,
                )
                and validate_delimited_file(
                    map_tmp,
                    min_columns=10,
                    min_data_rows=1,
                    has_header=False,
                )
            ):
                raise RuntimeError(
                    f"Segmentation outputs failed validation for {sample.name}"
                )

        update_stage_manifest(
            analysis_dict["output_dir"],
            f"segmentation:{sample.name}",
            segment_outputs,
            metadata=segment_metadata,
            input_paths=[sample.ratio_file],
        )

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
            posterior_list = []
            for item in merging_items:
                if int(item["start"]) < min_start:
                    min_start = int(item["start"])
                if int(item["end"]) > max_end:
                    max_end = int(item["end"])
                ratio_list.append(float(item["log2_ratio"]))
                region_list.append(item["region"])
                posterior_list.append(item["hmm_posterior"])

            mean_posterior = np.mean(posterior_list)
            mean_ratio = round(np.median(ratio_list), 3)
            new_segment = {
                "chr": first_dict["chr"],
                "start": min_start,
                "end": max_end,
                "region": ",".join(region_list),
                "n_regions": str(len(region_list)),
                "log2_ratio": mean_ratio,
                "state": first_dict["state"],
                "hmm_posterior": mean_posterior,
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
        posterior_list = []

        for item in merging_items:
            if int(item["start"]) < min_start:
                min_start = int(item["start"])
            if int(item["end"]) > max_end:
                max_end = int(item["end"])

            ratio_list.append(float(item["log2_ratio"]))
            region_list.append(item["region"])
            posterior_list.append(item["hmm_posterior"])

        mean_posterior = np.mean(posterior_list)
        mean_ratio = round(np.median(ratio_list), 3)
        new_segment = {
            "chr": first_dict["chr"],
            "start": min_start,
            "end": max_end,
            "region": ",".join(region_list),
            "n_regions": str(len(region_list)),
            "log2_ratio": mean_ratio,
            "state": first_dict["state"],
            "hmm_posterior": mean_posterior,
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
            coordinate = f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\t{tmp[3]}\t{tmp[4]}\t{tmp[5]}\t{tmp[-1]}\n"
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
                outline = f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\t{tmp[3]}\t{tmp[-1]}\n"
                o.write(outline)
            # f.close()
        o.close()

        rscript = sample.ratio_file.replace(".ratios.bed", ".CBS.R")
        segment_file = sample.ratio_file.replace(".bed", ".segment.bed")
        sample.add("segment_file", segment_file)

        r = open(rscript, "w")
        r.write("library(DNAcopy)" + "\n")
        line = f'cn <- read.table("{to_segment}", header=F)'
        r.write(line + "\n")
        line = "CNA.object <-CNA( genomdat = cn[,5], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')"
        r.write(line + "\n")
        line = "CNA.smoothed <- smooth.CNA(CNA.object)"
        r.write(line + "\n")
        line = f"segs <- segment(CNA.object, verbose=0, min.width={n_segments}, alpha = {alpha})"
        r.write(line + "\n")
        line = "segs2=segs$output"
        r.write(line + "\n")
        line = f'write.table(segs2[,2:6], file="{segment_file}",row.names=F, col.names=F, quote=F, sep="\t")'
        r.write(line + "\n")
        r.close()

        if not os.path.isfile(segment_file):

            msg = f" INFO: Segmenting sample {sample.name}"
            logging.info(msg)

            cmd = f"Rscript {rscript}"
            subprocess.run(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )

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
