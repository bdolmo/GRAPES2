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
    sample_object = ""
    for s in sample_list:
        if s.name != sample:
            name_tag = f"{s.name}_ratio"
            controls.append(name_tag)
        else:
            output_dir = s.sample_folder
            sample_object = s

    gene_plot_name = f"{sample}.{gene}.png"
    gene_plot = str(Path(output_dir) / gene_plot_name)

    df = pd.read_csv(analysis_dict["all_ratios"], sep="\t")
    df = df[df["exon"].str.endswith(gene)]
    df_dict = df.to_dict(orient="records")
    sample_tag = f"{sample}_ratio"

    controls_ratios = []
    exons_controls = []
    sample_ratios = []
    sample_exons = []

    all_ratios = []
    for idx,row in enumerate(df_dict):
        sample_exons.append(row["exon"]+";"+str(row["start"]))
        sample_ratios.append(float(row[sample_tag]))
        all_ratios.append(float(row[sample_tag]))
        for control in controls:
            if not control in row:
                continue
            controls_ratios.append(float(row[control]))
            exons_controls.append(row["exon"]+";"+str(row["start"]))
            all_ratios.append(float(row[control]))
    plot_dict = {"controls": controls_ratios, "exons": exons_controls}
    if not "per_exon" in sample_object.analysis_json:
        sample_object.analysis_json["per_exon"] = defaultdict(dict)

    sample_dict = dict()
    sample_dict[sample_tag] = sample_ratios
    sample_dict["exons"] = sample_exons
    sample_df = pd.DataFrame.from_dict(sample_dict)
    max_ratio = max(all_ratios) + 0.1

    min_ratio = -3.5

    # if min_ratio > -1:
    #     min_ratio = -1
    with sns.axes_style("dark"):
        fig, axes = plt.subplots(1, figsize=(22, 7))

        # DataFrame for boxplot
        plot_df = pd.DataFrame.from_dict(plot_dict)
        sample_df = pd.DataFrame.from_dict(sample_dict)

        # Merge and sort DataFrame
        combined_df = plot_df.merge(sample_df, on='exons')
        combined_df.sort_values(by='exons', inplace=True)

        sample_object.analysis_json["per_exon"][gene] = combined_df.to_json()

        # Plotting
        sns.set(font_scale=1.5)
        sns.set_style({'axes.linewidth': 0.5})
        fig, ax = plt.subplots(1, figsize=(22, 7))

        # Boxplot
        sns.boxplot(x="exons", y="controls", data=combined_df, showfliers=False, color="#d4ebf2", ax=ax)

        # Determine unique exons for x-axis po+-*/r4de xcsitions
        unique_exons = combined_df['exons'].unique()

        # Scatter plot for red and black points
        red_points = combined_df[combined_df[sample_tag] < analysis_dict["upper_del_cutoff"]]
        black_points = combined_df[(combined_df[sample_tag] >= analysis_dict["upper_del_cutoff"]) & (combined_df[sample_tag] <= analysis_dict["lower_dup_cutoff"])]

        # Plot each point separately to align with the corresponding box
        for exon in unique_exons:
            # Red points
            if exon in red_points['exons'].values:
                y_value = red_points[red_points['exons'] == exon][sample_tag].values[0]
                ax.scatter(unique_exons.tolist().index(exon), y_value, color='red', s=140)
            
            # Black points
            if exon in black_points['exons'].values:
                y_value = black_points[black_points['exons'] == exon][sample_tag].values[0]
                ax.scatter(unique_exons.tolist().index(exon), y_value, color='black', s=140)

        ax.axhline(0.433, color="blue")
        ax.axhline(-0.621, color="red")
        ax.set_xticklabels(unique_exons, rotation=90)
        ax.set(ylim=(-3.5, max_ratio))
        ax.set(ylabel="log2 ratio")

        for tick in ax.get_xticklabels():
            tick.set_color('black')
        for tick in ax.get_yticklabels():
            tick.set_color('black')

        fig.savefig(gene_plot, bbox_inches="tight")
        plt.close()

        return gene_plot, sample_object
        

def plot_single_exon_cnv(df, sample, variant_title):
    """Plots single exon CNVs with advanced visualization and aesthetics."""
    # Set the style
    sns.set(style="whitegrid", font_scale=1.2)
    
    # Create a color palette that's visually appealing and accessible
    primary_color = "#D55E00"  # A distinct color for the primary sample
    secondary_color_palette = sns.light_palette("#0072B2", len(df.columns[4:]), reverse=True)
    
    # Adjusting the figure size and resolution
    fig, ax1 = plt.subplots(figsize=(12, 7))

    # Plotting logic
    sample_cols = df.columns[4:]
    df["Position"] = df.index + df["start"].astype(int)
    
    for i, col in enumerate(sample_cols):
        if i == 0:  # Highlight the primary sample
            sns.lineplot(x=df["Position"], y=df[col], ax=ax1, label=col, color=primary_color, linewidth=2.5)
        else:  # Plot secondary samples with progressively lighter shades
            sns.lineplot(x=df["Position"], y=df[col], ax=ax1, label='_nolegend_', color=secondary_color_palette[i], linewidth=1.5, linestyle='--')
    
    ax1.set_ylabel("Log2 Ratio", fontsize=14)
    ax1.set_xlabel("Position", fontsize=14)
    ax1.set_title(variant_title, fontsize=16, fontweight='bold')
    ax1.set_ylim(-3.1, 1.25)

    # Enhance readability with custom ticks and gridlines
    ax1.grid(True, which='major', linestyle='--', linewidth='0.5', color='gray')
    ax1.legend(title="Controls", title_fontsize='13', fontsize='12', loc='upper right')
    
    # Secondary Y-axis for copy number
    ax2 = ax1.twinx()
    ax2.set_ylim(-3.1, 1.25)
    ax2.set_yticks([-3, -1, 0, 0.584, 1])
    ax2.set_yticklabels([0, 1, 2, 3, 4], fontsize=12)
    ax2.set_ylabel("Copy Number", fontsize=14)
    
    # Further style enhancements
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    
    sns.despine(right=False, ax=ax1)  # Removes the right spine for the primary y-axis
    sns.despine(ax=ax2, left=True)  # Removes the left spine for the secondary y-axis, keeping the right spine visible

    # Save the plot
    png_file = os.path.join(sample.sample_folder, f"{variant_title.replace(' ', '_')}.png")
    plt.savefig(png_file, format='png', dpi=300, bbox_inches='tight')

    # Clear the figure to free memory
    plt.close(fig)


def plot_svd_removal(df, png_file):

    plt.figure(figsize=(10, 10))
    sns.scatterplot(x='PC1', y='PC2', data=df)
    plt.title('First two Principal Components')
    plt.savefig(png_file, format='png', dpi=300)
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

    def plot_genomewide(self, genomewide, by_chr, sample):

        # Naming the genome plot
        plot_name = "{}{}".format(self._sample, ".genomewide.png")
        plot = str(Path(self._output_dir) / plot_name)


        cnr_df = pd.read_csv(self._cnr_file, sep="\t")
        cnr_df['chr'] = cnr_df['chr'].astype(str)  # Making sure the 'chr' column is string
        cnr_df = cnr_df.reindex(index=natsorted(cnr_df.index, key=lambda i: (cnr_df['chr'][i], cnr_df['start'][i])))
        cnr_df = cnr_df.reset_index(drop=True)  # Reset the index to maintain the new order

        sample_ratio = self._sample + "_ratio"
        min_ratio = cnr_df[sample_ratio].min()
        min_limit = -3.5
        if min_ratio <= -3:
            min_limit = -3.5

        max_ratio = cnr_df[sample_ratio].max()
        max_limit = 1
        if max_ratio > max_limit:
            max_limit = max_ratio+0.2
        
        # Setting chromosome color
        palette_dict = defaultdict(dict)
        color_list = ["#4f6b76", "#b2bbc0"]

        unique_chromosomes = natsorted(cnr_df["chr"].unique().tolist())
        idx = 0
        for chr in unique_chromosomes:
            if idx == 2:
                idx = 0
            palette_dict[chr] = color_list[idx]
            idx += 1
        x_list = []

        # Setting xtick divisions and labels
        chromosomes = natsorted(cnr_df["chr"].tolist())
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
            # if not os.path.isfile(plot):
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
            ratio_plot.set(ylim=(min_limit, max_limit))
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

            sample.analysis_json["genome_wide"] = cnr_df.to_json()

        return plot, sample
