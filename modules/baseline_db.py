import os
import sys
from sqlalchemy import create_engine, Column, Float, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import hashlib
import pandas as pd
import re

Base = declarative_base()

class BaseLine(Base):
    __tablename__ = 'BASELINE'
    id = Column(Integer, primary_key=True)
    sample = Column(String)
    bed = Column(String)
    bed_md5 = Column(String)
    chr = Column(String)
    start = Column(Integer)
    end = Column(Integer)
    exon = Column(String)
    gc = Column(Float)
    map = Column(Integer)
    normalized_depth = Column(Float)


def init_db(db_location):
    """
    Initialize the database and create tables.
    """
    db_location = f"sqlite:///{db_location}"
    engine = create_engine(db_location)
    Base.metadata.create_all(engine)
    return engine

def import_baselines_to_df(analysis_dict, df):
    """ """
    db_location = f'sqlite:///{analysis_dict["baseline_db"]}'
    engine = create_engine(db_location)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    if not analysis_dict["bed_md5"]:
        msg = f' ERROR: missing md5 for {analysis_dict["bed"]}'
        print(msg)

    baselines = session.query(BaseLine).filter_by(bed_md5=analysis_dict["bed_md5"]).all()
    
    # retrieve latest baseline number
    # unique_baselines = set()
    baselines_depth_dict = {}
    max_n = 0
    for baseline in baselines:
        result = re.search(r"\d+", baseline.sample)
        if result:
            baseline_n = int(result.group(0))
            if (baseline_n) > max_n:
                max_n = baseline_n
        if not baseline.sample in baselines_depth_dict:
            baselines_depth_dict[baseline.sample] = []
        baselines_depth_dict[baseline.sample].append(baseline.normalized_depth)
        # unique_baselines.add(baseline.sample)
    analysis_dict["latest_baseline"] = max_n

    for baseline_id in baselines_depth_dict:
        df[baseline_id] = baselines_depth_dict[baseline_id]

    return df


def calculate_baseline_median_depth(analysis_dict, sample_list, baselines):
    """
    Calculate the median _normalized_final depth for samples within the same baseline.
    """

    df = pd.read_csv(analysis_dict["normalized_depth"], sep='\t')
    # Ensure there's a column in df to match samples with baselines; assuming 'sample' here  
    # Add a column for median_normalized_final if it doesn't exist
   
    # Iterate through each baseline
    ref_idx = int(analysis_dict["latest_baseline"])+1

    if not baselines:
        for sample in sample_list:
            normalized_final_cols = [col for col in df.columns if '_normalized_final' in col]
            for col in normalized_final_cols:
                if col.startswith("baseline"):
                    df = df.drop(columns=[col])
                    continue
                new_ref_id = f"baseline{str(ref_idx)}_normalized_final"
                df[new_ref_id] = df[col]
                ref_idx+=1
        update_baselines(df, analysis_dict)
        return

    for baseline in baselines:
        # Filter the DataFrame to only include samples in the current baseline
        edit_baseline = ["chr","start","end","exon","gc","map"]
        for b in baseline:
            # edit_baseline.append(b)
            edit_baseline.append(f"{b}_normalized_final")
        baseline = edit_baseline
        baseline_df = df[baseline]

        baselines_n_samples = len(baseline)
        
        # Identify columns that contain _normalized_final values
        tmp_cols = [col for col in df.columns if '_normalized_final' in col]
        normalized_final_cols = []
        for col in tmp_cols:
            if col in baseline:
                normalized_final_cols.append(col)

        # Calculate the median of these _normalized_final values for each row
        median_values = baseline_df[normalized_final_cols].median(axis=1)

        new_ref_id = f"baseline{str(ref_idx)}_normalized_final"
        df[new_ref_id] = median_values
        for b in baseline:
            if b.startswith("baseline"):
                try:
                    df[new_ref_id]
                except:
                    pass
                else:
                    corr = df[b].corr(df[new_ref_id])
                    if corr >= 0.9999:
                        df = df.drop(columns=[new_ref_id])
                        continue
                df = df.drop(columns=[b])
            # edit_baseline.append(b)
        ref_idx+=1
        #sys.exit()
    update_baselines(df, analysis_dict)

def update_baselines(df, analysis_dict):
    """ """

    db_location = f'sqlite:///{analysis_dict["baseline_db"]}'
    engine = create_engine(db_location)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # get baseline referesces
    baseline_cols = [col for col in df.columns if 'baseline' in col]
 
    for baseline in baseline_cols:
        existing_entry = session.query(BaseLine).filter_by(sample=baseline).first()
        if not existing_entry:
            for index, row in df.iterrows():
                # Insert a new entry if it does not exist
                new_entry = BaseLine(
                    sample=baseline,  # Adjust accordingly
                    bed=analysis_dict["bed"],  # Adjust accordingly
                    bed_md5=analysis_dict["bed_md5"],  # Adjust accordingly or calculate
                    chr=row['chr'],
                    start=row['start'],
                    end=row['end'],
                    exon=row['exon'],
                    gc=row['gc'],
                    map=row['map'],
                    normalized_depth=row[baseline]
                )
                session.add(new_entry)
    session.commit()
    

def calculate_bed_md5(bed: str)-> str:
    """ """
    # create tmp bed
    tmp_bed = bed.replace(".bed", ".tmp.bed")
    o = open(tmp_bed, "w")
    with open(bed) as f:
        for line in f:
            tmp = line.split("\t")
            o.write(f"{tmp[0]}\t{tmp[1]}\t{tmp[2]}\n")
    o.close()

    md5_hash = hashlib.md5()
    # Open the file in binary mode
    with open(tmp_bed, 'rb') as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    os.remove(tmp_bed)

    return md5_hash.hexdigest()

