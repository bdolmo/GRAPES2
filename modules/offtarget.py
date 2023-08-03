import os
import pybedtools
import pandas as pd
from modules.utils import sort_bed_file
from modules.gc_content import annotate_gc_bed
from modules.mappability import annotate_mappability_bed







def extract_offtarget(sample, bed, genome_fasta, analysis_dict):
    """ """

    return sample



def create_offtarget_bed(bed, output_dir, genome_fasta, mappability_bed, chromosomes_file, blacklist):

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