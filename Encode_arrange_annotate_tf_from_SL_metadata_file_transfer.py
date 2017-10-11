#!/usr/bin/env python
import os
import re
import pandas as pd
import numpy as np
import glob

xls_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets.xlsx")
dir_path = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/SL*")
#dir_path = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/new_ones/SL*")
chip_analysis_path = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_motif_analysis_total/IDR*")
#chip_analysis_path = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_motif_analysis_total/new_ones/IDR*")

idr_file_list = glob.glob(dir_path)
idr_dir_list = glob.glob(chip_analysis_path)


### Read the excel file
#read_xls =  pd.read_excel(xls_file, sheetname = "Sheet1")  # or, sheetname = 0
read_xls = pd.ExcelFile(xls_file).parse("unique_TFs")
df_xls =  read_xls.iloc[:,[0,1,2,3,4]]
"""select the subset of xls rows that you want to operate on"""
#df_xls = df_xls.iloc[196:208]
#df_xls.shape[0]

df_xls.columns = ["rep1", "control1", "tf_name", "rep2", "control2",]
df_xls_ordered = df_xls.loc[:,["rep1","rep2","control1","control2","tf_name"]]
df_xls_ordered["tf_name"] = df_xls_ordered["tf_name"].apply(lambda x: "_".join(str(x).split("-")))
df_xls_ordered.to_csv("/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets.bed", header=True, sep="\t", index=True)

df_xls_ordered["combined_cols"] = df_xls_ordered.dropna().apply(lambda x: '_'.join(x.astype(str)), axis=1)
SL_all_split_list = df_xls_ordered[u'combined_cols'].dropna().str.split("_", 4).tolist()


SL_cross_check_list = [[],[]] # contains rep1_rep2 and tfname type list:
for each in SL_all_split_list:
    SL_cross_check_list[0].append("_".join(each[0:2])) # appends the combined the rep1 and rep2 :'SL88833_SL88834'
    SL_cross_check_list[1].append(each[4]) # appends the TF name


### Changes the name of files based on SLs i.e appends TFs:
for each_idr in idr_file_list:
    each_idr_basename= os.path.basename(each_idr)
    test_str = each_idr_basename   #test_str = 'SL101267.filt.nodup.srt.SE_VS_SL115566.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz'
    idr_string = re.split("[.,_]", test_str)  #idr_string = re.split("[.,_]", test_str)
    idr_SL_list = [ idr_string[i] for i in [0,6]] # Picks up the 0th SL and 6th SL element
    idr_check_string= "_".join(idr_SL_list)


    # compare the SL string from the idr_files to the SL no. from the excel file:
    for i,item in enumerate(SL_cross_check_list[0]):
        if idr_check_string == item :
            new_file_name =  each_idr + "_" + SL_cross_check_list[1][i]
            os.environ["old_file_name"]= each_idr
            os.environ["new_file_name"]= new_file_name
            os.system('mv $old_file_name $new_file_name')


### Changes the name of directory based on SLs, i.e appends TFs:
for each_idr in idr_dir_list:
    each_idr_basename= os.path.basename(each_idr)
    test_str = each_idr_basename   #test_str = 'SL101267.filt.nodup.srt.SE_VS_SL115566.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz'
    idr_string = re.split("[_]", test_str)  #idr_string = re.split("[.,_]", test_str)
    idr_SL_list = [ idr_string[i] for i in [1,2]]
    idr_check_string= "_".join(idr_SL_list)


    # compare the SL string from the idr_files to the SL no. from the excel file:
    for i,item in enumerate(SL_cross_check_list[0]):
        if idr_check_string == item:
            new_dir_name = each_idr + "_" + SL_cross_check_list[1][i]
            os.environ["old_dir_name"]= each_idr
            os.environ["new_dir_name"]= new_dir_name
            os.system('mv $old_dir_name $new_dir_name')


# find unannotated list of IDR dirs in chip_motif_analysis_total
# ls -I "IDR_*_*_*" 

# find unannotated list of IDR passed peak files in idr_passed_peaks_total
# ls *.gz






