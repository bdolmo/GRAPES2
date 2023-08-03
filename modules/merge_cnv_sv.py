import pandas as pd
import os
from pybedtools import BedTool
import re
def merge_bed_files(cnv_file, sv_file, output_file):

    """ """
    svtypes = ["DEL", "DUP"]

    # Check if files exist
    cnv_exists = os.path.isfile(cnv_file) and os.path.getsize(cnv_file) > 0
    sv_exists = os.path.isfile(sv_file) and os.path.getsize(sv_file) > 0

    if not cnv_exists and not sv_exists:
        raise FileNotFoundError("Both input files are missing.")
    elif cnv_exists and not sv_exists:
        # If only cnv_file exists, write it to output
        cnv_bed = pd.read_csv(cnv_file, sep='\t', names=["chr", "start", "end", "info"])
        cnv_bed = sort_bed_file(cnv_bed)
        cnv_bed.drop_duplicates(subset=["chr", "start", "end"], keep='first', inplace=True)
        cnv_bed.to_csv(output_file, sep='\t', header=False, index=False)
        return
    elif sv_exists and not cnv_exists:
        # If only sv_file exists, write it to output
        sv_bed = pd.read_csv(sv_file, sep='\t', names=["chr", "start", "end", "info"])
        sv_bed = sort_bed_file(sv_bed)
        sv_bed.drop_duplicates(subset=["chr", "start", "end"], keep='first', inplace=True)
        sv_bed.to_csv(output_file, sep='\t', header=False, index=False)
        return

    # If both files exist, proceed with the merging

    # Read the BED files using pandas
    cnv_bed = pd.read_csv(cnv_file, sep='\t', names=["chr", "start", "end", "info"])
    sv_bed = pd.read_csv(sv_file, sep='\t', names=["chr", "start", "end", "info"])

    # Create BedTool objects
    cnv_bt = BedTool.from_dataframe(cnv_bed)
    sv_bt = BedTool.from_dataframe(sv_bed)

    # Intersect BedTool objects
    intersected = cnv_bt.intersect(sv_bt, wa=True, wb=True)

    if intersected.count() == 0:
        cnv_bed = sort_bed_file(cnv_bed)
        sv_bed = sort_bed_file(sv_bed)
        combined_df = pd.concat([cnv_bed, sv_bed], ignore_index=True)
        combined_df.drop_duplicates(subset=["chr", "start", "end"], keep='first', inplace=True)
        combined_df.to_csv(output_file, sep='\t', header=False, index=False)
        return

    # Create a DataFrame from the intersected BedTool object
    intersected_df = intersected.to_dataframe(names=["chr", "start", "end", "cnv_info", "chr_sv", "start_sv", "end_sv", "sv_info"])
    # sys.exit()
    merged_df = pd.DataFrame(columns=intersected_df.columns)
    merged_df = sort_bed_file(merged_df)


    fields = ["SVTYPE", "REGION", "NREGIONS", "LOG2RATIO", "CN", "SCORE",
        "MAPQ", "KDIV", "BREAKREADS", "ASSEMBLED", "PE"]


    for svtype in svtypes:
        # Check that the SVTYPE is the same in both the CNV and SV files
        svtype_str = "SVTYPE=" + svtype
        temp_df = intersected_df[intersected_df.cnv_info.str.contains(svtype_str) & intersected_df.sv_info.str.contains(svtype_str)]
        
        # Use SV coordinates for the merged file
        temp_df['chr'] = temp_df['chr_sv']
        temp_df['start'] = temp_df['start_sv']
        temp_df['end'] = temp_df['end_sv']
        # Create new column 'info'
        temp_df['info'] = "UNDEF"

        # Merge the info fields, eliminating duplicates
        for idx, row in temp_df.iterrows():
            info_dict = {}
            for info_str in [row['cnv_info'], row['sv_info']]:
                info_items = info_str.split(';')
                # Loop over the items
                for item in info_items:
                    # Split the item into a key and value at the equals sign
                    if '=' in item:
                        key, value = item.split('=')
                        for field in fields:
                            if key == field:
                                info_dict[key] = value

                    # info_dict[key] = value

            temp_df.at[idx, 'info'] = "PRECISE" + ";"  +  ";".join([f"{k}={v}" for k, v in info_dict.items()])
        # Concatenate dataframes
        merged_df = pd.concat([merged_df, temp_df[['chr', 'start', 'end', 'info']]], ignore_index=True)



    # Write the merged BED file
    merged_df[["chr", "start", "end", "info"]].to_csv(output_file, sep='\t', header=False, index=False)


def sort_bed_file(df):
    if df.empty:
        return df
    # Replace 'chrX' and 'chrY' with temporary placeholders
    df['chr'] = df['chr'].replace({'chrX': 'chr23', 'chrY': 'chr24'})

    # Remove 'chr' prefix for sorting but keep it in a separate column
    df['chr_num'] = df['chr'].str.replace('chr', '').astype(int)

    # Sort by chromosomal position
    df = df.sort_values(['chr_num', 'start', 'end'])

    # Drop the temporary column
    df = df.drop(columns=['chr_num'])

    # Replace temporary placeholders with 'chrX' and 'chrY'
    df['chr'] = df['chr'].replace({'chr23': 'chrX', 'chr24': 'chrY'})

    return df
