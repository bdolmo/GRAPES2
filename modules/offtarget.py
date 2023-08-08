import os
import pybedtools
import pandas as pd
from modules.utils import sort_bed_file
from modules.gc_content import annotate_gc_bed
from modules.mappability import annotate_mappability_bed
import subprocess
from natsort import natsorted

def join_bins(input_bed, sample_names):


    sample_names = natsorted(sample_names)
    output_bed = input_bed.replace(".bed", ".tmp.bed")
    df = pd.read_csv(input_bed, sep="\t", header=None)

    num_samples = df.shape[1] - 6  # calculate number of sample columns
    sample_cols = [f"sample{i}" for i in range(1, num_samples + 1)]

    df.columns = ["chr", "start", "end", "exon", "gc", "map"] + sample_names

    df[['window', 'piece']] = df['exon'].str.split(';', expand=True)

    df['length'] = df['end'] - df['start']
    df['gc'] = df['gc'] * df['length']
    df['map'] = df['map'] * df['length']
    df = df[df['chr'] != 'chrM']  # remove rows where 'chrom' is 'chrM'
    # groupby and calculate mean for each sample
    grouped = df.groupby('window')

    result = pd.DataFrame({
        "chr": grouped['chr'].first(),
        "start": grouped['start'].min(),
        "end": grouped['end'].max(),
        "exon": grouped['exon'].first(),
        "gc": (grouped['gc'].sum() / grouped['length'].sum()).round(3),
        "map": (grouped['map'].sum() / grouped['length'].sum()).round(3)
    })

    for sample in sample_names:
        print(sample)
        result[sample] = grouped[sample].sum()

    result = result.reindex(natsorted(result.index))

    result.to_csv(output_bed, sep="\t", index=False, header=True)
    return output_bed


def extract_offtarget(sample_list, genome_fasta, ngs_utils_dict, analysis_dict):
    """ """
    
    msg = " INFO: Extracting off-target read depth information"
    print(msg)

    bams = natsorted([sample.bam for sample in sample_list ])
    sample_names = natsorted([sample.name for sample in sample_list])   

    offtarget_tmp_counts_name = "GRAPES2.offtarget.raw.counts.bed"
    offtarget_tmp_counts = f'{analysis_dict["output_dir"]}/{offtarget_tmp_counts_name}'

    # create a BedTool object from the bed file
    bed = pybedtools.BedTool(analysis_dict["offtarget_gc_map"])

    if not os.path.isfile(offtarget_tmp_counts):
        # calculate coverage
        result = bed.multi_bam_coverage(bams=bams)

        # save results to the output file
        result.saveas(offtarget_tmp_counts)

        tmp_bed = join_bins(offtarget_tmp_counts, sample_names)
        os.remove(offtarget_tmp_counts)
        os.rename(tmp_bed, offtarget_tmp_counts)

    analysis_dict["offtarget_raw_counts"] = offtarget_tmp_counts

    return analysis_dict



def create_offtarget_bed(bed, output_dir, genome_fasta, analysis_dict, mappability_bed, chromosomes_file, blacklist):

    # Load bed file into a BedTool object
    bedtool = pybedtools.BedTool(bed)

    # Define function to pad interval
    def pad_interval(interval):
        if interval.start < interval.stop:
            interval.start -= 1000
            interval.stop += 1000
        else:
            interval.start += 1000
            interval.stop -= 1000
        return interval

    # Apply padding to intervals
    bedtool = bedtool.each(pad_interval)

    # Save padded bed file
    padded_bed_path = f"{output_dir}/ontarget_paded_400bp.bed"
    bedtool.saveas(padded_bed_path)

    # Concatenate blacklist file to the padded bed file
    concat_bed_path = f"{output_dir}/ontarget_paded_400bp_centromers_patches.bed"
    bedtool = pybedtools.BedTool(padded_bed_path)\
        .cat(pybedtools.BedTool(blacklist), postmerge=False)
    
    # Exclude intervals with "_" in their name
    bedtool = bedtool.filter(lambda interval: '_' not in interval.chrom)

    # Save concatenated bed file
    bedtool.saveas(concat_bed_path)
    sort_bed_file(concat_bed_path)

    bedtool = pybedtools.BedTool(concat_bed_path)

    # Complementary file generation for the peak detection step
    complement_bedtool = bedtool.complement(g=chromosomes_file)
    
    # Filter intervals by length and add "Complement" as name
    complement_bedtool = complement_bedtool.filter(lambda interval: interval.length >= 1)
    complement_bedtool = complement_bedtool.each(lambda interval: pybedtools.create_interval_from_list([interval.chrom, str(interval.start), str(interval.end), "Complement"]))
    
    # Save off-target bed file
    complement_bedtool.saveas(f"{output_dir}/offtarget_unprocessed.bed")
    sort_bed_file(f"{output_dir}/offtarget_unprocessed.bed")

    # Pseudo offtargets
    offtarget_windows_bed = create_pseudowindows(output_dir, 1000000, chromosomes_file)

    gc_bed = annotate_gc_bed(offtarget_windows_bed, output_dir, genome_fasta)
    map_bed = annotate_mappability_bed(gc_bed, mappability_bed)

    analysis_dict["offtarget_gc_map"] = map_bed

    return analysis_dict
    # Additiona annotations such as GC and mappability


def create_pseudowindows(output_dir, binsize, chromosomes_file):
    """ """
    chromosomes_dict = {}
    with open(chromosomes_file) as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            chr = tmp[0]
            end = int(tmp[1])
            chromosomes_dict[chr] = end
    f.close()


    inputpw  = os.path.join(output_dir, 
        "offtarget_unprocessed.bed")

    outputpw = os.path.join(output_dir, 
        f"offTarget.pseudowindowed.bed")

    try:
        INPW = open(inputpw, 'r')
    except IOError:
        print(f"Cannot open {inputpw}")
        return

    try:
        OUTPW = open(outputpw, 'w')
    except IOError:
        print(f"Cannot open {outputpw}")
        return

    sum     =  0
    presum  =  0
    lc      =  0
    winc    =  1
    prechr  = ''
    cntpiece=  1
    newend  = ''
    newst   = ''
    prend   = ''
    prest   = ''

    for line in INPW:
        line = line.rstrip()
        lc += 1
        tmp = line.split('\t')
        chr = tmp[0]
        st  = int(tmp[1])
        end = int(tmp[2])

        if st == 0 and end == 0:
            continue

        if lc==1: 
            prechr=chr

        length = end-st
        sum += length

        if chr != prechr:
            prechr=chr
            winc += 1
            cntpiece = 1
            sum = length
            presum = 0

        if presum == sum and sum == binsize:
            prend = st
            presum = 0
            continue

        if sum > binsize:
            diff = binsize-presum
            iters = round(sum/binsize)
            for i in range(iters):
                if i == 0:
                    newst = st
                    newend = newst + diff
                    sum = binsize

                    if newst >= chromosomes_dict.get(chr, 0):
                        continue

                    if newend >= chromosomes_dict.get(chr, 0):
                        newend = chromosomes_dict.get(chr, 0)

                    OUTPW.write(f"{chr}\t{newst}\t{newend}\tpdwindow_{winc};piece_{cntpiece}\n")

                    cntpiece = 1
                    winc += 1
                    prend = newend
                    prest = newst
                    presum = sum

                    nlen = end - prend
                    if nlen <= binsize:
                        sum = nlen
                        newst = prend + 1
                        newend = end

                        if newend > chromosomes_dict.get(chr, 0):
                            newend = chromosomes_dict.get(chr, 0)

                        OUTPW.write(f"{chr}\t{newst}\t{newend}\tpdwindow_{winc};piece_{cntpiece}\n")

                        cntpiece += 1
                        presum = sum
                        break
                    else:
                        sum = binsize
                        newst = prend + 1
                        newend = newst + binsize

                        if newend > chromosomes_dict.get(chr, 0):
                            newend = chromosomes_dict.get(chr, 0)

                        OUTPW.write(f"{chr}\t{newst}\t{newend}\tpdwindow_{winc};piece_{cntpiece}\n")

                        cntpiece = 1
                        winc += 1
                        prend = newend
                        prest = newst
                        continue
                else:
                    nlen = end - prend
                    if nlen <= binsize:
                        sum = nlen
                        newst = prend + 1
                        newend = end
                        if newst >= chromosomes_dict.get(chr, 0):
                            continue
                        if newend > chromosomes_dict.get(chr, 0):
                            newend = chromosomes_dict.get(chr, 0)

                        OUTPW.write(f"{chr}\t{newst}\t{newend}\tpdwindow_{winc};piece_{cntpiece}\n")

                        cntpiece += 1
                        prend = newend
                        prest = newst
                        presum = sum
                        break
                    else:
                        newst = prend + 1
                        newend = newst + binsize
                        sum = binsize
                        if newst >= chromosomes_dict.get(chr, 0):
                            continue
                        if newend > chromosomes_dict.get(chr, 0):
                            newend = chromosomes_dict.get(chr, 0)

                        OUTPW.write(f"{chr}\t{newst}\t{newend}\tpdwindow_{winc};piece_{cntpiece}\n")

                        cntpiece = 1
                        winc += 1
                        prend = newend
                        prest = newst
                        continue
            continue

        if end > chromosomes_dict.get(chr, 0):
            end = chromosomes_dict.get(chr, 0)

        OUTPW.write(f"{chr}\t{st}\t{end}\tpdwindow_{winc};piece_{cntpiece}\n")

        if sum == binsize:
            cntpiece = 1
            winc += 1
            prend = end
            prest = st
            presum = sum
            continue
        else:
            cntpiece += 1
            prend = end
            prest = st
            presum = sum
            continue

        presum = sum
        if presum == binsize:
            presum = 0

    INPW.close()
    OUTPW.close()

    return outputpw