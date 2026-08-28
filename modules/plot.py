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


def _genomewide_scatter_style(point_count):
    """Return a readable marker size and opacity for the target density."""
    if point_count >= 100000:
        return 3.0, 0.18
    if point_count >= 50000:
        return 4.0, 0.25
    if point_count >= 10000:
        return 7.0, 0.40
    return 14.0, 0.70


def _genomewide_y_limits(ratios, del_cutoff, dup_cutoff):
    """Choose robust limits without allowing rare artifacts to flatten the plot."""
    finite_ratios = pd.Series(ratios)[np.isfinite(ratios)]
    if finite_ratios.empty:
        return -1.5, 1.5

    lower_quantile = finite_ratios.quantile(0.001)
    upper_quantile = finite_ratios.quantile(0.999)
    lower_limit = min(-1.5, float(del_cutoff) - 0.4, lower_quantile - 0.15)
    upper_limit = max(1.5, float(dup_cutoff) + 0.4, upper_quantile + 0.15)

    # Values beyond these bounds are still shown as triangles at the plot edge.
    return max(-3.5, lower_limit), min(3.5, upper_limit)


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
        cnr_df[sample_ratio] = pd.to_numeric(cnr_df[sample_ratio], errors="coerce")
        finite_mask = np.isfinite(cnr_df[sample_ratio])
        finite_ratios = cnr_df.loc[finite_mask, sample_ratio]
        if finite_ratios.empty:
            msg = (
                " WARNING: No finite ratio values found for sample %s while "
                "plotting genomewide profile; using default y-axis limits"
            )
            logging.warning(msg, self._sample)
        min_limit, max_limit = _genomewide_y_limits(
            finite_ratios,
            self._del_cutoff,
            self._dup_cutoff,
        )
        
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
        # Put chromosome labels at their centers and separators at boundaries.
        chromosome_starts = []
        chromosome_centers = []
        for chromosome in unique_chromosomes:
            chromosome_indices = cnr_df.index[cnr_df["chr"] == chromosome]
            chromosome_starts.append(int(chromosome_indices[0]))
            chromosome_centers.append(
                (int(chromosome_indices[0]) + int(chromosome_indices[-1])) / 2
            )

        if genomewide == True:
            sns.set_style("ticks")
            fig, ratio_plot = plt.subplots(figsize=(20, 7), dpi=180)
            ratio_plot.set_title(self._sample, fontsize=20, pad=16)

            finite_x = cnr_df.index.to_numpy()[finite_mask]
            finite_y = finite_ratios.to_numpy()
            finite_colors = (
                cnr_df.loc[finite_mask, "chr"].map(palette_dict).to_numpy()
            )
            marker_size, marker_alpha = _genomewide_scatter_style(len(finite_y))
            lower_outliers = finite_y < min_limit
            upper_outliers = finite_y > max_limit
            visible_points = ~(lower_outliers | upper_outliers)

            ratio_plot.scatter(
                finite_x[visible_points],
                finite_y[visible_points],
                c=finite_colors[visible_points],
                s=marker_size,
                alpha=marker_alpha,
                edgecolors="none",
                linewidths=0,
                rasterized=True,
            )

            clipped_marker_size = max(8.0, marker_size * 1.8)
            if upper_outliers.any():
                ratio_plot.scatter(
                    finite_x[upper_outliers],
                    np.full(upper_outliers.sum(), max_limit - 0.03),
                    c=finite_colors[upper_outliers],
                    marker="^",
                    s=clipped_marker_size,
                    alpha=max(0.55, marker_alpha),
                    edgecolors="none",
                    linewidths=0,
                    rasterized=True,
                )
            if lower_outliers.any():
                ratio_plot.scatter(
                    finite_x[lower_outliers],
                    np.full(lower_outliers.sum(), min_limit + 0.03),
                    c=finite_colors[lower_outliers],
                    marker="v",
                    s=clipped_marker_size,
                    alpha=max(0.55, marker_alpha),
                    edgecolors="none",
                    linewidths=0,
                    rasterized=True,
                )

            # Setting y limits
            ratio_plot.set(ylim=(min_limit, max_limit))
            ratio_plot.set_xlim(-0.5, max(len(cnr_df) - 0.5, 0.5))
            ratio_plot.set_xticks(chromosome_centers)
            ratio_plot.set_xticklabels(
                unique_chromosomes,
                rotation=45,
                rotation_mode="anchor",
                ha="right",
                size=12,
            )
            ratio_plot.tick_params(axis="y", labelsize=15)
            ratio_plot.yaxis.set_major_locator(ticker.MaxNLocator(nbins=8))
            ratio_plot.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
            ratio_plot.set_xlabel("", fontsize=20)
            ratio_plot.set_ylabel("log2 Ratio", fontsize=20)
            ratio_plot.axhline(0, color="#687078", lw=0.7, alpha=0.55, zorder=0)
            ratio_plot.axhline(
                self._dup_cutoff,
                ls="--",
                color="#1769e0",
                lw=1.1,
                alpha=0.9,
            )
            ratio_plot.axhline(
                self._del_cutoff,
                ls="--",
                color="#d62728",
                lw=1.1,
                alpha=0.9,
            )
            ratio_plot.grid(axis="y", color="#d7dce0", lw=0.6, alpha=0.55)
            ratio_plot.grid(axis="x", visible=False)

            # Adding vertical lines to separate chromosomes
            for chromosome_start in chromosome_starts[1:]:
                ratio_plot.axvline(
                    x=chromosome_start - 0.5,
                    color="#525b63",
                    lw=0.6,
                    alpha=0.45,
                    zorder=0,
                )

            clipped_count = int(lower_outliers.sum() + upper_outliers.sum())
            if clipped_count:
                ratio_plot.text(
                    0.995,
                    0.985,
                    (
                        f"Triangles: {int(upper_outliers.sum())} above / "
                        f"{int(lower_outliers.sum())} below display range"
                    ),
                    transform=ratio_plot.transAxes,
                    ha="right",
                    va="top",
                    fontsize=9,
                    color="#4d545a",
                )

            # Saving as png
            fig.tight_layout()
            fig.savefig(plot, dpi=180, bbox_inches="tight", facecolor="white")
            plt.close(fig)

            sample.analysis_json["genome_wide"] = cnr_df.to_json()

        return plot, sample
