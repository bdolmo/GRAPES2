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
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from adjustText import adjust_text

pd.options.mode.chained_assignment = None  # default='warn'


def plot_normalization(sample_list, analysis_dict):
    """
    Plotting coverage normalization
    """

    df = pd.read_csv(analysis_dict["normalized_depth"], sep="\t")

    for sample in sample_list:

        msg = (" INFO: Plotting normalized coverage for sample {}").format(sample.name)
        logging.info(msg)

        png_name = ("{}.normalization.png").format(sample.name)
        normalization_plot = str(Path(analysis_dict["output_dir"]) / png_name)
        sample.add("normalization_plot", normalization_plot)
        if not os.path.isfile(normalization_plot):

            normalized_gc_tag = ("{}_normalized_gc").format(sample.name)
            normalized_map_tag = ("{}_normalized_map").format(sample.name)

            sns.set(font_scale=2)
            fig, axes = plt.subplots(2, 2, figsize=(25, 22))
            fig.suptitle("GC-content & Mappability correction", fontsize=50)
            axes[0, 0] = sns.boxplot(
                ax=axes[0, 0],
                x="gc_integer",
                y=sample.name,
                data=df,
                showfliers=False,
                palette="Blues",
            )
            axes[0, 0].set_title("Raw coverage vs GC", fontsize=30)
            axes[0, 0].set_xticklabels(axes[0, 0].get_xticklabels(), rotation=30)
            axes[0, 0].set(xlabel="%GC")
            axes[0, 0].set(ylabel="Coverage")
            axes[0, 0].xaxis.set_major_locator(ticker.MultipleLocator(base=5))

            axes[0, 1] = sns.boxplot(
                ax=axes[0, 1],
                x="gc_integer",
                y=normalized_gc_tag,
                data=df,
                showfliers=False,
                palette="Blues",
            )
            axes[0, 1].set_title("GC-content corrected coverage", fontsize=30)
            axes[0, 1].set_xticklabels(axes[0, 1].get_xticklabels(), rotation=30)
            axes[0, 1].set(xlabel="%GC")
            axes[0, 1].set(ylabel="Coverage")
            axes[0, 1].xaxis.set_major_locator(ticker.MultipleLocator(base=5))

            axes[1, 0] = sns.boxplot(
                ax=axes[1, 0],
                x="map_integer",
                y=sample.name,
                data=df,
                showfliers=False,
                palette="Blues",
            )
            axes[1, 0].set_title("Raw coverage vs mappability", fontsize=30)
            axes[1, 0].set_xticklabels(axes[1, 0].get_xticklabels(), rotation=30)
            axes[1, 0].set(xlabel="%Mappability")
            axes[1, 0].set(ylabel="Coverage")
            axes[1, 0].set(xlim=(0, 100))
            axes[1, 0].xaxis.set_major_locator(ticker.MultipleLocator(base=10))

            axes[1, 1] = sns.boxplot(
                ax=axes[1, 1],
                x="map_integer",
                y=normalized_map_tag,
                data=df,
                showfliers=False,
                palette="Blues",
            )
            axes[1, 1].set_title("GC-Mappability corrected coverage", fontsize=30)
            axes[1, 1].set_xticklabels(axes[1, 1].get_xticklabels(), rotation=30)
            axes[1, 1].set(xlabel="%Mappability")
            axes[1, 1].set(ylabel="Coverage")
            axes[1, 1].set(xlim=(0, 100))
            axes[1, 1].xaxis.set_major_locator(ticker.MultipleLocator(base=10))
            fig.savefig(normalization_plot)

    return sample_list


def plot_gene(sample, sample_list, gene, analysis_dict):
    """ """

    controls = []
    output_dir = ""
    for s in sample_list:
        if s.name != sample:
            name_tag = ("{}_ratio").format(s.name)
            controls.append(name_tag)
        else:
            output_dir = s.sample_folder

    gene_plot_name = ("{}.{}.png").format(sample, gene)
    gene_plot = str(Path(output_dir) / gene_plot_name)
    if os.path.isfile(gene_plot):
        return
    df = pd.read_csv(analysis_dict["all_ratios"], sep="\t")
    df = df[df["exon"].str.endswith(gene)]
    df_dict = df.to_dict(orient="records")

    sample_tag = ("{}_ratio").format(sample)

    controls_ratios = []
    exons_controls = []
    sample_ratios = []
    sample_exons = []

    for row in df_dict:
        sample_exons.append(row["exon"])
        sample_ratios.append(float(row[sample_tag]))

        for control in controls:
            controls_ratios.append(float(row[control]))
            exons_controls.append(row["exon"])
    plot_dict = {"controls": controls_ratios, "exons": exons_controls}
    sample_dict = dict()
    sample_dict[sample_tag] = sample_ratios
    sample_dict["exons"] = sample_exons
    sample_df = pd.DataFrame.from_dict(sample_dict)
    max_ratio = max(sample_ratios) + 0.1
    min_ratio = min(sample_ratios) - 0.1

    if max_ratio < 0.5:
        max_ratio = 0.5
    if min_ratio > -1:
        min_ratio = -1
    with sns.axes_style("dark"):
        fig, axes = plt.subplots(1, figsize=(22, 7))
        plot_df = pd.DataFrame.from_dict(plot_dict)
        sns.set(font_scale=1.5)

        ratio_tag = ("{}_ratio").format(sample)

        axes = sns.boxplot(
            x="exons", y="controls", data=plot_df, showfliers=False, color="#d4ebf2"
        ).set_title(gene)
        axes = sns.scatterplot(
            x="exons", y=sample_tag, data=sample_df, s=140, color="red"
        )
        axes.axhline(0.433, ls="--", color="blue")
        axes.axhline(-0.621, ls="--", color="red")
        axes.set_xticklabels(axes.get_xticklabels(), rotation=90)
        axes.set(ylim=(min_ratio, max_ratio))
        axes.set(ylabel="log2 ratio")

        fig.savefig(gene_plot, bbox_inches="tight")
        plt.close()


def plot_single_exon(data, case_sample, exon, output_file):
    """ """
    df = pd.DataFrame.from_dict(data)

    min_ratio = df["log2_ratio"].min()
    max_ratio = df["log2_ratio"].max()

    min_limit = -1
    if min_ratio < min_limit:
        min_limit = min_ratio - 0.2

    max_limit = 0.58
    if max_ratio > max_limit:
        max_limit = max_limit + 0.2

    palette = {
        c: "red" if c == case_sample else "darkgrey" for c in df["sample"].unique()
    }
    sns.set(rc={"figure.dpi": 180, "savefig.dpi": 180})
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(15, 7))
    # fig.suptitle(self._sample, fontsize=20)
    ax.axhline(0.433, ls="--", color="blue", zorder=1)
    ax.axhline(-0.621, ls="--", color="red", zorder=1)

    ax = sns.lineplot(
        data=df,
        x="position",
        y="log2_ratio",
        hue="sample",
        size="class",
        linewidth=1.5,
        sizes=(3, 1.5),
        palette=palette,
        legend=False,
        zorder=2,
    )

    label_format = "{:,.0f}"
    ticks_loc = ax.get_xticks().tolist()
    ax.set_xticks(ax.get_xticks().tolist())
    ax.set_xticklabels([label_format.format(x) for x in ticks_loc])

    # after plotting the data, format the labels
    # current_values = ax.get_xticks()
    # using format string '{:.0f}' here but you can choose others
    # ax.set_xticklabels(['{:.0f}'.format(x) for x in current_values])

    # ax.set_xticklabels(ax.get_xticklabels(),rotation =60)
    ax.set_title(case_sample + "-" + exon, fontsize=20)
    ax.set(ylim=(min_limit, max_limit))
    ax.set_ylabel("log2 Ratio", fontsize=15)
    ax.set_xlabel("coordinate", fontsize=15)

    # ax.set_xticks(ticks, labels, rotation=45, ha='right', rotation_mode='anchor')
    fig.savefig(output_file)
    plt.close()


class CnvPlot:
    """
    Class for plotting cnvs
    """

    def __init__(
        self, cnr_file, cns_file, calls, sample, output_dir, dup_cutoff, del_cutoff
    ):
        # ratios
        self._cnr_file = cnr_file
        # segments
        self._cns_file = cns_file
        # calls
        self._calls = calls
        # sample name
        self._sample = sample
        # output directory
        self._output_dir = output_dir

        self._dup_cutoff = dup_cutoff

        self._del_cutoff = del_cutoff

    def plot_genomewide(self, genomewide, by_chr):

        # Naming the genome plot
        plot_name = "{}{}".format(self._sample, ".genomewide.png")
        plot = str(Path(self._output_dir) / plot_name)

        if os.path.isfile(plot):
            return plot

        cnr_df = pd.read_csv(self._cnr_file, sep="\t")

        sample_ratio = self._sample + "_ratio"
        min_ratio = cnr_df[sample_ratio].min()
        min_limit = -1
        if min_ratio < min_limit:
            min_limit = min_ratio - 0.2

        # Setting chromosome color
        palette_dict = defaultdict(dict)
        color_list = ["#4f6b76", "#b2bbc0"]

        unique_chromosomes = cnr_df["chr"].unique().tolist()
        idx = 0
        for chr in unique_chromosomes:
            if idx == 2:
                idx = 0
            palette_dict[chr] = color_list[idx]
            idx += 1
        x_list = []

        # Setting xtick divisions and labels
        chromosomes = cnr_df["chr"].tolist()
        i = 0
        xtick_list = []
        for chr in chromosomes:
            if not chr in x_list:
                x_list.append(chr)
                xtick_list.append(i)
            else:
                x_list.append("")
            i += 1

        if genomewide == True:
            if not os.path.isfile(plot):
                sns.set(rc={"figure.dpi": 180, "savefig.dpi": 180})
                sns.set_style("ticks")
                fig, axes = plt.subplots(figsize=(20, 7))
                fig.suptitle(self._sample, fontsize=20)

                ratio_plot = sns.scatterplot(
                    data=cnr_df,
                    x=cnr_df.index,
                    y=cnr_df[sample_ratio],
                    size=0.1,
                    hue=cnr_df.chr,
                    palette=palette_dict,
                    edgecolor="none",
                )

                # Setting y limits
                ratio_plot.set(ylim=(min_limit, 1))
                ratio_plot.set_xticks(xtick_list)
                ratio_plot.set_xticklabels(unique_chromosomes, rotation=40, size=15)
                ratio_plot.set_yticks(ratio_plot.get_yticks())
                ratio_plot.set_yticklabels(ratio_plot.get_yticks(), size=15)
                ratio_plot.get_legend().remove()
                ratio_plot.set_xlabel("", fontsize=20)
                ratio_plot.set_ylabel("log2 Ratio", fontsize=20)
                ratio_plot.axhline(self._dup_cutoff, ls="--", color="blue")
                ratio_plot.axhline(self._del_cutoff, ls="--", color="red")

                # Adding vertical lines to separate chromosomes
                for xc in xtick_list:
                    ratio_plot.axvline(x=xc, color="black")
                # Saving as png
                ratio_plot.figure.savefig(plot)
                plt.close()
        return plot
