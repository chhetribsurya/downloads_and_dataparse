import pandas as pd
import pybedtools
import os
import re
from os.path import join
from os.path import basename
from os.path import splitext


output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/metadata_hepg2_dl"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

snyder_file = "/gpfs/gpfs1/home/schhetri/for_SLs/snyders_final_fastq_download_metadata.txt"
brads_file = "/gpfs/gpfs1/home/schhetri/for_SLs/brads_final_fastq_download_metadata.txt"
stam_peggy_file = "/gpfs/gpfs1/home/schhetri/for_SLs/stam_peggy_final_fastq_download_metadata.txt"

snyder_df = pd.read_csv(snyder_file, sep="\t")
snyder_df = snyder_df.drop("Unnamed: 0", axis=1)
snyder_df["lab"] = "Snyder"
snyder_df

brads_df = pd.read_csv(brads_file, sep="\t")
brads_df = brads_df.drop("Unnamed: 0", axis=1)
brads_df["lab"] = "Bradstein"
brads_df

stam_peggy_df = pd.read_csv(stam_peggy_file, sep="\t")
stam_peggy_df = stam_peggy_df.drop("Lab_1", axis=1)
stam_peggy_df = stam_peggy_df.drop("Unnamed: 0", axis=1)
stam_peggy_df = stam_peggy_df.rename(columns={"Lab_2":"lab"})
stam_peggy_df

# comibine the dataframes of all labs:
combined_df = pd.concat([snyder_df, brads_df, stam_peggy_df], ignore_index=True)
combined_tf_df = combined_df[~combined_df["Experiment_target_1"].str.contains(r"^H2|H3|H4")].reset_index().drop("index", axis=1)
combined_histone_df = combined_df[combined_df["Experiment_target_1"].str.contains(r"^H2|H3|H4")].reset_index().drop("index", axis=1)

combined_tf_df.to_excel(join(output_dir, "Encode_full_hepg2_datasets.xls"), sheet_name="Sheet1", header=True, index=True)
combined_tf_df.to_csv(join(output_dir, "Encode_full_hepg2_datasets.txt"), sep="\t", header=True, index=True)

combined_histone_df.to_excel(join(output_dir, "Encode_hepg2_chipseq_histones_dl_data.xls"), sheet_name="Sheet1", header=True, index=True)
combined_histone_df.to_csv(join(output_dir, "Encode_hepg2_chipseq_histones_dl_data.txt"), sep="\t", header=True, index=True)

