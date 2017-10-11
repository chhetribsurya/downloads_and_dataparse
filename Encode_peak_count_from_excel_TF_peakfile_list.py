#!/usr/bin/env python

import pandas as pd
import glob
import os
import subprocess
import re
from os.path import join
from os.path import basename
from collections import defaultdict


#dir_path = os.path.expanduser("~/for_chris/batch_I/idr_passed_peaks_total/SL*")
#file_list = glob.glob(dir_path)

output_dir = os.path.expanduser("~/for_chris/hepg2_xls_info")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

excel_file = os.path.expanduser("~/for_chris/Encode_full_hepg2_datasets.xlsx")
# df_xls = pd.read_excel(excel_file, sep="\t", sheet_name="Sheet2")


"""Read multiple excel sheets from the same file"""
xls = pd.ExcelFile(excel_file)
df_xls = xls.parse("Sheet1")
# df_xls2 = xls.parse("Main_sheet"); df_xls3 = xls.parse("dups"); df_xls4 = xls.parse("flagged")

df_xls["Target"] = df_xls["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
TF_list = df_xls["Target"]

### cutting down the df to make it manageable:
df_xls_select =df_xls.iloc[:,0:6]


### Calculates the peak count of excel file row by rows maintaining the order:
file_path = "~/for_chris/batch_I/idr_passed_peaks_total"
peak_count_list = []
for index,row in df_xls_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name

	full_path = os.path.expanduser(join(file_path, SL_rep1 + "*" + SL_rep2 + "*"))
	file_name = glob.glob(full_path)[0]
	print file_name
	num_lines = sum(1 for line in open(file_name))
	print tf_name, num_lines, lab_name, "\n"
	peak_count_list.append(num_lines)

### maintains the same order as the original excel file:
df_xls_select["peak_count"] = peak_count_list
df_xls_select.to_csv(join(output_dir, "Encode_full_hepg2_datasets_with_peakcount"), header=True, sep="\t", index=True)
df_xls_select.to_excel(join(output_dir, "Encode_full_hepg2_datasets_with_peakcount.xls"), header=True, sheet_name="Sheet1")


### ordered or sorted by the TF name and lab name of your interest ( lab order custom selected):
concat_list = []
df_myers = df_xls_select[df_xls_select["Lab"].str.contains("Myers")]
df_myers = df_myers.sort_values("Target")
concat_list.append(df_myers)

df_snyder = df_xls_select[df_xls_select["Lab"].str.contains("Snyder")]
df_snyder = df_snyder.sort_values("Target")
concat_list.append(df_snyder)

df_bradstein = df_xls_select[df_xls_select["Lab"].str.contains("Bradstein")]
df_bradstein = df_bradstein.sort_values("Target")
concat_list.append(df_bradstein)

df_peggy = df_xls_select[df_xls_select["Lab"].str.contains("Peggy ")]
df_peggy = df_peggy.sort_values("Target")
concat_list.append(df_peggy)

df_Iyer = df_xls_select[df_xls_select["Lab"].str.contains(" Iyer")]
df_Iyer = df_Iyer.sort_values("Target")
concat_list.append(df_Iyer)

combined_df = pd.concat(concat_list)
combined_df.to_csv(join(output_dir, "Encode_full_hepg2_datasets_with_peakcount_labsorted"), header=True, sep="\t", index=True)
combined_df.to_excel(join(output_dir, "Encode_full_hepg2_datasets_with_peakcount_labsorted.xls"), header=True, sheet_name="Sheet1")



#######################################
#######################################
#######################################
#######################################


"""Find all the duplicated lines based on TF names, and print the no. of peak counts"""

### make a default dict to find the duplicates 
### based on item name( i.e TF name here):
tf_dup_dict = defaultdict(list)
for i,item in enumerate(TF_list):
	tf_dup_dict[item].append(i)


### For bunch printing of duplicate lines:
dup_idx = []
for key,values in tf_dup_dict.iteritems():
	if len(values) >1:
		print key, ":", values
		dup_idx.extend(values)

df_xls_dups = df_xls_select.iloc[dup_idx]


### For one by one assessment of the duplicates within, and
### do certain operations for each of duplicated idx values:
dup_peak_count_list = []
for key,values in tf_dup_dict.iteritems():
	if len(values) > 1:
		for each_idx in values:
			#print df_xls_select.iloc[each_idx]
			SL_rep1 = df_xls_select.iloc[each_idx]["Exp ID rep1"]
			SL_rep2 = df_xls_select.iloc[each_idx]["Exp ID rep2"]
			tf_name = df_xls_select.iloc[each_idx]["Target"]
			lab_name = df_xls_select.iloc[each_idx]["Lab"]
			print SL_rep1, SL_rep2, tf_name, lab_name

			file_path = "~/for_chris/batch_I/idr_passed_peaks_total"
			full_path = os.path.expanduser(join(file_path, SL_rep1 + "*" + SL_rep2 + "*"))
			file_name = glob.glob(full_path)[0]
			print file_name
			num_lines = sum(1 for line in open(file_name))
			print tf_name, num_lines, lab_name, "\n"
			dup_peak_count_list.append(num_lines)

df_xls_dups["peak_count"] = dup_peak_count_list
df_xls_dups
df_xls_dups.to_csv(join(output_dir, "Encode_full_hepg2_datasets_dups"), header=True, sep="\t", index=True)
df_xls_dups.to_excel(join(output_dir, "Encode_full_hepg2_datasets_dups.xls"), header=True, sheet_name="Sheet1")

# general environment variables, preparing for pipeline run:   
# os.environ["peak_file_full_path"] = each_tf_file_in_file_list_from_glob
# output_peakcount = subprocess.check_output('wc -l $peak_file_full_path | cut -f1 -d " " ', shell=True)
# print tf_name, ":", output_peakcount






