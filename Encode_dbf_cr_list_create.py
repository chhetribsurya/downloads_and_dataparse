import pandas as pd
import numpy as np
import pybedtools
import pickle
import scipy 
import os
import re
import glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from rpy2.robjects.packages import importr  #stats = importr('stats')
from rpy2.robjects.vectors import FloatVector
from os.path import join
from os.path import basename
from os.path import splitext


main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/tf_co_binding_total"
if not os.path.exists(main_dir):
	os.makedirs(main_dir)

""" Choose the file category or TF category you want to analyse on """
#file_category = "all_tf" #file_category = "all_tf"; file_category = "cr_tf"; file_category = "dbf_tf"

""" All TF containing dir """
#dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
all_tf_file_list = glob.glob(dir_path) # represents all tf file list

###################################

"""Read multiple excel sheets from the same file"""
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets_DBF_CR.xls"
df_xls = pd.ExcelFile(excelfile).parse("unique_TFs_DBF_CR_annotation")
df_xls["Category_split_join"] = df_xls["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
xls_tf_list =  df_xls["Target"].tolist()

cr_df = df_xls[df_xls["Category"] == "CR/CF"]
dbf_df = df_xls[~(df_xls["Category"] == "CR/CF")]

""" Check the xls TF list and the file suffix tf_list to make sure, if they are in the TF list"""
TF_check_list = []
for each_file in all_tf_file_list:
	tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
	TF_check_list.append(tf_name)

# df_xls[df_xls["Target"].isin(TF_check_list)] # or,
if sorted(TF_check_list) == sorted(xls_tf_list):
	print "\nGreat!! TF list resembles b/w the file list and xls tf list..."
else:
	print "\nWarning: Your files TF list doesn't resemble the xls tf list..."


"""Select the category of TF for chosing the file list"""
"""Gives the path list containing the file path for each TF"""
cr_tf_file_list = []
for each_tf in cr_df["Target"].tolist():
	for each_file in all_tf_file_list:
		if each_file.endswith(each_tf):
			cr_tf_file_list.append(each_file)

dbf_tf_file_list = []
for each_tf in dbf_df["Target"].tolist():
	for each_file in all_tf_file_list:
		if each_file.endswith(each_tf):
			dbf_tf_file_list.append(each_file)


print "DBF and CR/CF list created!! You could iterate on each TF for desired operation"



