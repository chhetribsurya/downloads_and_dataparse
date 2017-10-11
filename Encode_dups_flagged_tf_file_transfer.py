#!/usr/bin/env python

import pandas as pd
import glob
import os
import subprocess
import shutil
from distutils.dir_util import copy_tree
import re
from os.path import join
from os.path import basename
from collections import defaultdict



output_dir = os.path.expanduser("~/for_chris/batch_I/idr_passed_peaks_total")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sub_dir_dup = join(output_dir, "dups")
if not os.path.exists(sub_dir_dup):
	os.makedirs(sub_dir_dup)


sub_dir_flag = join(output_dir, "flagged")
if not os.path.exists(sub_dir_flag):
	os.makedirs(sub_dir_flag)


sub_dir_unique = join(output_dir, "unique_TFs")
if not os.path.exists(sub_dir_unique):
	os.makedirs(sub_dir_unique)


meme_output_dir = os.path.expanduser("~/for_chris/batch_I/chip_motif_analysis_total")
if not os.path.exists(meme_output_dir):
    os.makedirs(meme_output_dir)

sub_meme_dir_flagged = join(meme_output_dir, "flagged")
if not os.path.exists(sub_meme_dir_flagged):
	os.makedirs(sub_meme_dir_flagged)

sub_meme_dir_dup = join(meme_output_dir, "dups")
if not os.path.exists(sub_meme_dir_dup):
	os.makedirs(sub_meme_dir_dup)

sub_meme_dir_unique = join(meme_output_dir, "unique_TFs")
if not os.path.exists(sub_meme_dir_unique):
	os.makedirs(sub_meme_dir_unique)

### The excel file which is processed with duplicates and flagged details:
### Output of script peak_count_of_peak_files_from_excel_TFlist.py
excel_file = os.path.expanduser("~/for_chris/Encode_full_hepg2_datasets.xlsx")

"""Read the duplicate file from multiple excel sheets of the same file"""
df_xls = pd.ExcelFile(excel_file).parse("dups")
df_xls["Target"] = df_xls["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
df_xls = df_xls.dropna()

### cutting down the df to make it manageable:
df_xls_select = df_xls.iloc[:,0:6]


""" Duplicate files and dir transfer"""

### Transfer the duplicate TF files using row by rows TF_name and SL#:
peak_file_path = "~/for_chris/batch_I/idr_passed_peaks_total"

for index,row in df_xls_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(peak_file_path, SL_rep1 + "*" + SL_rep2 + "*" ))
	file_name = glob.glob(full_path)[0]; #print file_name
	shutil.move(file_name, sub_dir_dup)


### Transfer the duplicate motif dir:
motif_file_path = "~/for_chris/batch_I/chip_motif_analysis_total"

for index,row in df_xls_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(motif_file_path, "IDR" + "*" + SL_rep1 + "*" + SL_rep2 + "*" ))
	dir_name = glob.glob(full_path)[0]; #print dir_name
	shutil.move(dir_name, sub_meme_dir_dup)




######################
######################
######################



"""Read the flagged file from multiple excel sheets of the same file"""
df_xls_flagged = pd.ExcelFile(excel_file).parse("flagged")
df_xls_flagged["Target"] = df_xls_flagged["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
df_xls_flagged = df_xls_flagged.dropna()

### cutting down the df to make it manageable:
df_xls_flagged_select = df_xls_flagged.iloc[:,0:6]

""" Flagged files and dir transfer"""

### Transfer the flagged TF files using row by rows TF_name and SL#:
peak_file_path = "~/for_chris/batch_I/idr_passed_peaks_total"

for index,row in df_xls_flagged_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(peak_file_path, SL_rep1 + "*" + SL_rep2 + "*" ))
	file_name = glob.glob(full_path)[0] ; #print file_name
	shutil.move(file_name, sub_dir_flag)


### Transfer the flagged motif dir:
motif_file_path = "~/for_chris/batch_I/chip_motif_analysis_total"

for index,row in df_xls_flagged_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(motif_file_path, "IDR" + "*" + SL_rep1 + "*" + SL_rep2 + "*" ))
	dir_name = glob.glob(full_path)[0]; #print dir_name
	shutil.move(dir_name, sub_meme_dir_flagged)



#########################
#########################
#########################

def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)


"""Read the unique_TF file from multiple excel sheets of the same file"""
df_xls = pd.ExcelFile(excel_file).parse("unique_TFs")
df_xls["Target"] = df_xls["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
df_xls = df_xls.dropna() # df_xls = df_xls.dropna(how="all")

### cutting down the df to make it manageable:
df_xls_unique = df_xls.iloc[:,0:6]


""" unique files and dir transfer"""

### Transfer the unique TF files using row by rows TF_name and SL#:
peak_file_path = "~/for_chris/batch_I/idr_passed_peaks_total"

for index,row in df_xls_unique.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(peak_file_path, SL_rep1 + "*" + SL_rep2 + "*" ))
	file_name = glob.glob(full_path)[0]; #print file_name
	shutil.copy(file_name, sub_dir_unique)


### Transfer the unique motif dir:
motif_file_path = "~/for_chris/batch_I/chip_motif_analysis_total"

for index,row in df_xls_unique.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(motif_file_path, "IDR" + "*" + SL_rep1 + "*" + SL_rep2 + "*" ))
	dir_name = glob.glob(full_path)[0]; #print dir_name
	
	os.environ["dir_name"] = dir_name
	os.environ["output_dir"] = sub_meme_dir_unique
	os.system("cp -r $dir_name $output_dir")

	



	
