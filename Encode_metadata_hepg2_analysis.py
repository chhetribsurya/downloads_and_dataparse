#!/bin/bash
import os
import re
import pandas as pd
import numpy as np
import glob

xls_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/HepG2_missing_SLs.xls.xlsx")
dir_path = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/idr_passed_peaks/SL*")

read_xls =  pd.read_excel(xls_file, sheetname = 0)
df_xls =  read_xls.iloc[:,[0,1,2,3,4,6]]

Target_col_list = [u'Exp ID rep1', u'Control ID rep1', u'Target', u'Exp ID rep2', u'Control ID rep2']
df_xls_select = df_xls.loc[:, df_xls.columns.isin(Target_col_list)]
df_xls_select.columns = ["rep1", "control1", "tf_name", "rep2", "control2",]
df_xls_ordered = df_xls_select.loc[:,["rep1","rep2","control1","control2","tf_name"]]

df_xls["combined_col"] = df_xls_ordered.dropna().apply(lambda x: '_'.join(x.astype(str).values), axis=1)
SL_all_split_list = df_xls[u'combined_col'].dropna().str.split("_", 4).tolist()

SL_cross_check_list = [[],[]]
for each in SL_all_split_list:
    SL_cross_check_list[0].append("_".join(each[0:2]))
    SL_cross_check_list[1].append(each[4])

#####################################################

idr_list = glob.glob(dir_path)

#tf_filepath_metadata = {}

count = 0
tf_filepath_metadata = [[],[],[]]
for each_idr in idr_list:
    each_idr_basename= os.path.basename(each_idr)
    test_str = each_idr_basename   #test_str = 'SL101267.filt.nodup.srt.SE_VS_SL115566.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz'
    #idr_string = filter(None, re.split("[.,_]", test_str))   #idr_string = re.split("[.,_]", test_str)
    idr_string = re.split("[.,_]", test_str)  #idr_string = re.split("[.,_]", test_str)
    idr_SL_list = [ idr_string[i] for i in [0,6]]
    idr_check_string= "_".join(idr_SL_list)

    for i,item in enumerate(SL_cross_check_list[0]):
        if idr_check_string == item :
            count += 1
            tf_filepath_metadata[0].append(SL_cross_check_list[1][i])
            tf_filepath_metadata[1].append(each_idr)
            tf_filepath_metadata[2].append(idr_check_string)
            print(count, SL_cross_check_list[1][i])

print ("\n\n")

my_analysis_list = df_xls["Prelim-Analysis on"].dropna().tolist()
for each in my_analysis_list:
    index = tf_filepath_metadata[0].index(each)
    print(tf_filepath_metadata[0][index],"\t",tf_filepath_metadata[1][index],"\t",tf_filepath_metadata[2][index])