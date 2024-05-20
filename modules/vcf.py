import pysam
from datetime import datetime
from natsort import natsorted, ns
import pybedtools
from pybedtools import BedTool
import os
import pysam


def annotate_snv_baf(bam, ref_fasta, chrom, start, end, cnv_type, copy_number) -> float:
    samfile = pysam.AlignmentFile(bam, "rb")
    fastafile = pysam.FastaFile(ref_fasta)  # Load the reference genome
    snv_allele_frequencies = []
    average_af = 0

    snv_dict = []
    total_compatible_vars = 0
    perc_compatible = 0
    # Use mpileup to analyze the region covered by the CNV/SV
    for pileupcolumn in samfile.pileup(chrom, start-1, end, stepper='samtools', min_mapping_quality=50, min_base_quality=20):
        if start <= pileupcolumn.pos + 1 <= end:
            pos = pileupcolumn.pos
            ref_base = fastafile.fetch(chrom, pos, pos+1).upper()
            total_count = 0
            variant_counts = {}
            has_indels_or_skips = False

            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    has_indels_or_skips = True
                    break 
                if pileupread.indel > 1:
                    has_indels_or_skips = True
                    break

            if has_indels_or_skips:
                continue

            for pileupread in pileupcolumn.pileups:
                # try:
                if pileupread.alignment.query_qualities[pileupread.query_position] > 20:
                    base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    if base != ref_base:  # Check if the base is different from the reference
                        variant_counts[base] = variant_counts.get(base, 0) + 1
                    total_count += 1
                # except:
                #     print("ERROR", bam, pileupread)
                #     sys.exit()
                #     pass

            # Calculate AF for each detected variant within the CNV/SV
            for base, count in variant_counts.items():
                if total_count > 0:
                    af = round(count / total_count, 3)

                    if total_count > 30 and count > 10:  # Apply filtering criteria for variant calling

                        if cnv_type == "DUP" and af == 1:
                            continue

                        is_compatible = is_vaf_compatible_with_cnv(af, cnv_type, copy_number)
                        var_dict = {
                            "chromosome": chrom,
                            "pos": pos,
                            "ref": ref_base,
                            "alt": base,
                            "depth": total_count,
                            "allele_depth": count,
                            "alele_frequency": af,
                            "compatible": is_compatible
                        }

                        if af >= 0.2:
                            if is_compatible:
                                total_compatible_vars+=1

                            snv_dict.append(var_dict)
                            snv_allele_frequencies.append(af)

    # Calculate and annotate average AF if any SNVs were found
    if snv_allele_frequencies:
        average_af = round(sum(snv_allele_frequencies) / len(snv_allele_frequencies),3)

        perc_compatible = round((total_compatible_vars/len(snv_dict)),3)

    # Close files
    samfile.close()
    fastafile.close()
    return snv_dict, average_af, perc_compatible


def is_vaf_compatible_with_cnv(vaf, cnv_type, copy_number=2):
    """
    Check if a Variant Allele Frequency (VAF) is compatible with a given Copy Number Variation (CNV),
    taking into account specific expectations for deletions and duplications.
    """
    tolerance = 0.125  # Adjust tolerance as necessary

    if cnv_type == 'DEL':
        # For a deletion, expect the VAF to be at least 0.8
        return vaf >= 0.8
    elif cnv_type == 'DUP':

        # Define expected VAFs based on the number of extra copies
        if copy_number == 3:  # One extra copy
            expected_vaf = 1 / 3
            if abs(vaf-expected_vaf) <= tolerance:
                return True
            expected_vaf = 2 / 3
            if abs(vaf-expected_vaf) <= tolerance:
                return True
            else:
                return False    
        elif copy_number == 4:  # Two extra copies
            expected_vaf = 2 / 4
        
        return abs(vaf - expected_vaf) <= tolerance



def get_reference_contigs_from_bam(bam_file_path):
    """
    Retrieve a list of reference contigs from the header of a BAM file.

    :param str bam_file_path: input bam file
    :returns: A list of strings, each representing a reference contig name.
    :rtype: list
    """
    full_contigs = []
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        contigs = [ref for ref in bam_file.references]
        header = dict(bam_file.header)
        reference_genome = header.get('SQ', [{}])[0].get('UR', 'Unknown')

        for contig in contigs:
            contig_length = bam_file.get_reference_length(contig)
            full_contigs.append(f"{contig},length={contig_length}")

    return full_contigs, reference_genome


def define_genotype_based_on_cn(cn: int):
    """ """

    if cn == 0:
        return '1/1'  # Homozygous deletion
    elif cn == 1:
        return '0/1'  # Heterozygous deletion
    elif cn == 2:
        return '0/0'  # Homozygous reference
    elif cn == 3:
        return '0/1'  # Heterozygous duplication
    elif cn >= 4:
        # For CN=4, it could be considered a homozygous duplication in some contexts
        # For CN > 4, the representation might depend on the interpretation and the specific variant caller
        return '1/1'  # Homozygous duplication or higher copy gain
    else:
        return '.'  # Undefined or invalid CN


def define_genotype_based_on_af(af: float):
    """ """
    if af >= 0.9:
        return '1/1'
    else:
        return '0/1'


def annotate_ontarget_overlaps(bed, roi_bed):
    """
    Tag each bed entry with an annotation (ON_TARGET=1 if it overlaps or ON_TARGET=0 if not) within the ROI bed and write to a file
    """
    a = pybedtools.BedTool(bed)
    b = pybedtools.BedTool(roi_bed)

    # Intersect and append the overlap width at the end of each record from 'a'
    result = a.intersect(b, wao=True)

    output_bed = bed.replace(".bed", ".tmp.bed")
    seen_records = {}
    with open(output_bed, 'w') as out_f:
        for feature in result:
            str_feature = str(feature)
            tmp_feature = str_feature.split("\t")

            output_line = "\t".join(tmp_feature[0:4])

            # feature[-1] contains the overlap size
            overlap_size = int(feature[-1])

            # Add 'ON_TARGET' field based on overlap
            annotated_feature = output_line
            if not "ON_TARGET" in annotated_feature:
                if overlap_size > 0:
                    annotated_feature = output_line + ";ON_TARGET=1\n"
                else:
                    annotated_feature = output_line + ";ON_TARGET=0\n"
            
            if not annotated_feature in seen_records:
                out_f.write(annotated_feature)
            seen_records[annotated_feature] = seen_records

    out_f.close()
    os.remove(bed)
    os.rename(output_bed, bed)
    

def bed_to_vcf(bed, roi_bed, bam, ref_fasta, output_vcf, sample, min_gc, max_gc, min_map, min_size):
    """
    Convert a BED file to a VCF file, using reference contigs extracted from the given BAM file to populate the VCF header.

    This function processes each line of the BED file, creating VCF entries with the structural variant information provided. 
    It writes a VCF file with appropriate headers and formatted entries.

    :param str bam: input bam file
    :param str bed: input bed with all merged CNV and SV calls
    :param str output_vcf: exported vcf
    :returns: None
    """
    now = datetime.now()
    date_time = now.strftime("%Y%m%d")

    reference_contigs, reference_genome = get_reference_contigs_from_bam(bam)
    vcf_header = [
        '##fileformat=VCFv4.3',
        '##source=GRAPES2',
        f'##reference={reference_genome}',
        f'##fileDate={date_time}'
    ]

    sample.analysis_json["calls"] = []

    # Add the reference genome contig information to the header
    for contig in reference_contigs:
        vcf_header.append(f'##contig=<ID={contig}>')

    # Add the column header
    info_fields = [
        '##FILTER=<ID=Incompatible_BAF,Description="Mean B-Allele Frequency is not compatible with the CNV">',
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Size of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">',
        '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
        '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
        '##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Mapping quality of SV">',
        '##INFO=<ID=KDIV,Number=1,Type=Float,Description="K-mer divergence">',
        '##INFO=<ID=REGION,Number=1,Type=String,Description="Name region of the 4th field in the BED file">',
        '##INFO=<ID=NREGIONS,Number=1,Type=Integer,Description="Total number of baits affected">',
        '##INFO=<ID=LOG2RATIO,Number=1,Type=Float,Description="Log2 CNV ratio">',
        '##INFO=<ID=GC,Number=1,Type=Float,Description="Mean GC-content of the affected region">',
        '##INFO=<ID=MAP,Number=1,Type=Float,Description="Mean mappability of the affected region">',
        '##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy Number">',
        '##INFO=<ID=BREAKREADS,Number=1,Type=Integer,Description="Number of reads supporting the breakpoint">',
        '##INFO=<ID=ASSEMBLED,Number=1,Type=Integer,Description="Number of reads assembled at breakpoints">',
        '##INFO=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the SV">',
        '##INFO=<ID=ZSCORE,Number=1,Type=Float,Description="Z-score of the CNV call">',
        '##INFO=<ID=CNV_SCORE,Number=1,Type=Float,Description="Phred score of the CNV call">',
        '##INFO=<ID=MBQ,Number=1,Type=Integer,Description="Mean Base Quality in phred scale">',
        '##INFO=<ID=MBAF,Number=1,Type=Float,Description="Mean B-Allele Frequency for overlapping SNV">',
        '##INFO=<ID=SNV,Number=1,Type=String,Description="Total overlapping SNVs">',
        '##INFO=<ID=BAF_COMPAT_SNV,Number=1,Type=String,Description="Ratio of BAF-compatible SNVs with an overlapping CNV">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">',
        '##INFO=<ID=CALLER,Number=1,Type=String,Description="Caller name">',
        '##INFO=<ID=ON_TARGET,Number=1,Type=Integer,Description="On-target call(1) or off-target(0)">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Variant Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    ]
    for field in info_fields:
        vcf_header.append(field)

    annotate_ontarget_overlaps(bed, roi_bed)

    bed_records = []
    bed_dict = {}

    # Read the BED file and process each line
    with open(bed, 'r') as bed_file:
        for line in bed_file:
            line = line.rstrip("\n")
            bed_records.append(line)      
    bed_file.close()

    bed_roi_records = []
    with open(roi_bed, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            tmp = line.split("\t")
            chrom = tmp[0]
            start = tmp[1]
            end = tmp[2]
            bed_roi_records.append((chrom, int(start), int(end)))
    f.close()

    tmp_input_bed = bed.replace(".bed", ".tmp.bed")
    newo = open(tmp_input_bed, "w")

    # Write the header to the VCF file
    bed_records = natsorted(bed_records)
    vcf_header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE')
    with open(output_vcf, 'w') as vcf_file:
        vcf_file.write('\n'.join(vcf_header) + '\n')
        for line in bed_records:

            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            end = fields[2]
            filter_tag = "."

            info_fields = fields[3].split(';')
            svtype = [f for f in info_fields if "SVTYPE=" in f][0].split('=')[1]
            genotype = "./."
            baf_compat_snv = "."
            cn = ""
            for field in info_fields:
                if field.startswith("CN="):
                    cn = int(field.replace("CN=", ""))
                    genotype = define_genotype_based_on_cn(cn)
                if field.startswith("AF="):
                    af = float(field.replace("AF=", ""))
                    genotype = define_genotype_based_on_af(af)

            snv_list = []
            baf = "."
            if cn:
                # annotate overlapping SNV
                snv_list, baf, baf_compat_snv = annotate_snv_baf(bam, ref_fasta, chrom, int(pos), int(end), svtype, cn)

                if cn == 1 and len(snv_list) > 0 and baf_compat_snv < 0.8:
                    filter_tag = "Incompatible_BAF"

            str_snv_list = []
            out_snv = "."
            for snv_dict in snv_list:
                str_snv_list.append(str(snv_dict))
            if str_snv_list:
                out_snv = ','.join(str_snv_list)

            ci_pos, ci_end = [0, 0], [0, 0]  # Default CI values
            # Check if the variant is within any BED region
            for bed_chrom, bed_start, bed_end in bed_roi_records:
                if chrom == bed_chrom and bed_start <= int(pos) <= bed_end:
                    # Calculate CIPOS and CIEND based on the variant position and region boundaries
                    ci_pos = [int(pos) - bed_start, int(pos) - bed_start]
                    ci_end = [int(pos) - bed_end, int(pos) - bed_end]
                    break
            if fields[3].startswith("PRECISE"):
                ci_pos = (-10, 10)
                ci_end = (-10, 10)
            if fields[3].startswith("IMPRECISE"):
                for idx,field in enumerate(info_fields):
                    if field.startswith("KDIV="):
                        info_fields[idx] = "KDIV=."
                    if field.startswith("MBQ="):
                        info_fields[idx] = "MBQ=."
            svlen = int(end)-int(pos)
            if svlen < min_size:
                continue
            
            # Construct the INFO field for VCF
            info_fields.insert(2, f"BAF_COMPAT_SNV={str(baf_compat_snv)}")
            info_fields.insert(2, f"SNV={str(out_snv)}")
            info_fields.insert(2, f"MBAF={str(baf)}")
            info_fields.insert(2, "CALLER=GRAPES2")
            info_fields.insert(2, "SOURCE=GRAPES2")
            info_fields.insert(2, f"CIEND=0,{ci_end[1]}")
            info_fields.insert(2, f"CIPOS={ci_pos[0]},0")
            info_fields.insert(2, f"SVLEN={svlen}")
            info_fields.insert(2, f"END={end}")

            call_dict = {
                "coordinates": f"{chrom}:{pos}-{end}",          
            }
            filter_by_gc = False
            filter_by_size = False
            filter_by_map = False

            for idx,field in enumerate(info_fields):
                tmp_field = field.split("=")
                field_name = tmp_field[0]
                
                if "PRECISE" in field:
                    call_dict["precision"] = field
                else:
                    if len(tmp_field) > 1:
                        field_value = tmp_field[1]
                        call_dict[field_name] = field_value.replace(";", "_")
                if field_name=="GC":
                    field_value = tmp_field[1]
                    if field_value != ".":
                        if float(field_value) < min_gc:
                            filter_by_gc = True
                        if float(field_value) > max_gc:
                            filter_by_gc = True
                info_fields[idx] = field.replace(";", "_")

            if filter_by_gc:
                continue
            if filter_by_map:
                continue

            call_dict["GENOTYPE"] = genotype
            sample.analysis_json["calls"].append(call_dict)

            info = ';'.join(info_fields)

            if filter_tag == ".":
                newo.write(f"{chrom}\t{pos}\t{end}\t{info}\n")


            # Create the VCF entry
            alt = f'<{svtype}>'
            vcf_entry = f'{chrom}\t{pos}\t.\tN\t{alt}\t.\t{filter_tag}\t{info}\tGT\t{genotype}\n'
            vcf_file.write(vcf_entry)
    vcf_file.close()
    newo.close()

    os.remove(bed)
    os.replace(tmp_input_bed, bed)

    return sample

