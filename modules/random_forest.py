import re
import pickle
import pandas as pd
from pathlib import Path

main_dir = Path(__file__).resolve().parents[1]
ml_models_dir = os.path.join(main_dir, "ml_models")

def extract_info_field(info):
    """
        Parses the INFO field and extracts relevant features.
    """
    info_dict = {}
    for entry in info.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[entry] = True  # Flags with no value
    
    # Convert numerical fields to float or int
    for key in ["GC", "MAP", "ZSCORE", "CV", "NREGIONS", "LOG2RATIO", "CNV_SCORE"]:
        if key in info_dict:
            try:
                info_dict[key] = float(info_dict[key])
            except ValueError:
                info_dict[key] = None
    
    return info_dict

def load_model(pickle_file="trained_rf.pkl"):
    """
        Loads the pre-trained RandomForest model from a pickle file.
    """

    pickle_file = os.path.join(ml_models_dir, "trained_rf.pkl")

    with open(pickle_file, "rb") as f:
        model = pickle.load(f)
    return model

def process_vcf(input_vcf, sample_correlation, sample_enrichment, output_vcf, model):
    """Reads the VCF file, applies the classifier, and saves the modified VCF."""

    # input_vcf = "/raw-data/projectes/GRAPES/GRAPES_ML/REGULOME_SEQ/GRAPES2/NA12249.nodup.simulated/NA12249.nodup.simulated.GRAPES2.vcf"
    # output_vcf = input_vcf.replace(".vcf", "rf.vcf")

    # feature_cols = ["gc", "map", "log2ratio", "cv", "nregions", "%ROI", "svlen", "mean_top10_corr"]
    feature_cols = [
        "svlen", "gc", "map", "zscore", "cv", "log2ratio", "nregions",
        "%ROI", "mean_top10_corr"
    ]

    updated_lines = []
    
    with open(input_vcf, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                updated_lines.append(line.strip())  # Preserve header lines
                continue

            fields = line.strip().split("\t")
            chrom, start, end, info = fields[0], fields[1], fields[2], fields[7]          

            # Extract features from INFO field
            info_dict = extract_info_field(info)

            end = info_dict["END"]

            info_dict_translated = {
                "svlen": int(end)-int(start),
                "gc": info_dict["GC"],
                "map": info_dict["MAP"],
                "log2ratio": info_dict["LOG2RATIO"],
                "cv": info_dict["CV"],
                "%ROI": sample_enrichment,
                "nregions": info_dict["NREGIONS"],
                "mean_top10_corr": sample_correlation,
                "zscore": info_dict["ZSCORE"]
            }

            # Prepare input data for prediction
            feature_values = [info_dict_translated.get(col, 0) for col in feature_cols] 

            feature_df = pd.DataFrame([feature_values], columns=feature_cols)

            # Apply model prediction
            prob_score = model.predict_proba(feature_df)[0, 1]  

            # Append prediction to INFO field
            new_info = f"{info};RF_SCORE={prob_score:.4f}"
            fields[7] = new_info  # Update INFO field
            updated_lines.append("\t".join(fields))

    # Save modified VCF
    with open(output_vcf, "w") as outfile:
        outfile.write("\n".join(updated_lines) + "\n")
    
    print(f"Updated VCF saved to: {output_vcf}")

if __name__ == "__main__":
    input_vcf = "input.vcf"
    output_vcf = "output_with_rf.vcf"

    # Load the trained model
    rf_model = load_model("trained_rf.pkl")

    # Process VCF and add RF_SCORE
    process_vcf(input_vcf, output_vcf, rf_model)
