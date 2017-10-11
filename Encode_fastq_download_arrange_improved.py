import pandas as pd
import pybedtools
import glob
import re
import os
import numpy
from itertools import izip
from os.path import join
from os.path import splitext
from os.path import basename

#metadata_input_file = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_snyder_fastq_metadata.tsv"
#metadata_input_file = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/downoad_bradstein_fastq_metadata.tsv" 
#metadata_input_file = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_stam_peggy_fastq_metadata.tsv"
metadata_input_file = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_IV_libraries/downoad_bradstein_fastq_metadata.tsv"

#output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_fastqs/snyder"
#output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_fastqs/bradstein"
#output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_fastqs/stam_peggy"
output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_IV_libraries/download_fastqs/brad"

#lab_name = "Snyder" #For all lab; lab_name = ""
#lab_name = "Bernstein" #For all lab; lab_name = ""
lab_name = "" #For all lab; lab_name = ""

batch_name = "batch_IV"
SL_start_num = 100 # Change this value, if you need unique SLs# for different batches as SL00001 could overlap b/w batches.


if not os.path.exists(output_dir):
	os.makedirs(output_dir)

read_file = pd.read_csv(metadata_input_file, sep="\t")
imp_cols = ["File accession", "Experiment accession", "Controlled by", "Biological replicate(s)", "Experiment target", "Biosample term name", \
			"Read length", "Lab", "File download URL",  "File Status",  "Audit NOT_COMPLIANT", "Audit INTERNAL_ACTION", 'Audit WARNING', "md5sum"]
read_file = read_file.loc[:, imp_cols]
select_cols = ["File accession", "Controlled by", "Experiment target", "Biosample term name", "Biological replicate(s)", "Read length", "Experiment accession", "Lab"]
read_df = read_file.loc[:, select_cols]
read_df.columns = read_df.columns.str.replace(" ", "_")

#control_list = read_df['Controlled by'].dropna().replace(r".*files\/(.*)\/", r"\1", regex=True).unique().tolist()
control_list = read_df['Controlled_by'].dropna().replace(r".*files\/([0-9A-Za-z]+)\/*", r"\1", regex=True).unique().tolist()
reps_list = read_df['File_accession'].unique().tolist()

#check if the fastq was already downloaded:
downloaded_file_list = [ basename(each_fastq) for each_fastq in glob.glob(join(output_dir,"*fastq.gz")) ]

### Download the controls if not available with in the fastq:
for each in control_list:
	if each not in reps_list:
		print "Control Fastq missing:", each+".fastq.gz"
		if each+".fastq.gz" not in downloaded_file_list:
			print "Downloading control fastq.....:", each
			os.environ["fastq_id"] = each
			os.environ["output_dir"] = output_dir
			os.system("wget --directory-prefix=${output_dir} https://www.encodeproject.org/files/${fastq_id}/@@download/${fastq_id}.fastq.gz")
		else:
			print "But, no problem, %s already downloaded...." %(each+".fastq.gz")




### Select the samples with only controls available:
read_df = read_df[~read_df['Controlled_by'].isnull()]
read_df["Controlled_by"] = read_df['Controlled_by'].replace(r".*files\/([0-9A-Za-z]+)\/*", r"\1", regex=True)

### Selecting the lab to be analysed on:
read_df = read_df[read_df['Lab'].str.lower().str.contains(lab_name.lower())]

### Select the samples with 2 or more replicates:
rep_idx = []
tf_name_list = read_df.loc[:,"Experiment_target"].tolist()
for i,item in enumerate(tf_name_list):
	if (tf_name_list.count(item) >= 2):
		print(item), i 
		rep_idx.append(i)

final_rep_df = read_df.iloc[rep_idx]
final_rep_df = final_rep_df.reset_index()
tf_list = final_rep_df.loc[:,"Experiment_target"].tolist()

### Edit all the tfs with more than 2 replicate counts with additional suffix:
new_tf_list = []
for item in tf_list:
	if new_tf_list.count(item) >= 2:
		new_tf_list.append(item + "_2")
	else:
		new_tf_list.append(item)

final_rep_df["Experiment_target"] = new_tf_list

### Pick all those tfs with 2 rep counts since some of those might have three rep counts, and after rename 
### those 3 set replicated could have lone replicates so, gotta eliminate those ones from the tf list:
### For instance: 'MXI1-human', 'MXI1-human', 'MXI1-human_2' (lone replicate)
final_rep_idx = []
for i,item in enumerate(new_tf_list):
	if (new_tf_list.count(item) != 1):
		print(item),i
		final_rep_idx.append(i)

final_df = final_rep_df.iloc[final_rep_idx]
final_read_df = final_df.loc[:,["File_accession", "Controlled_by", "Experiment_target"]]


rep1_df = final_read_df.iloc[::2]
rep2_df = final_read_df.iloc[1::2]

combined_df = pd.concat([rep1_df.reset_index(),rep2_df.reset_index()], axis=1)


if rep1_df.shape == rep2_df.shape:
	print "\n\nRep1 and Rep2 dataframe are of equal sizes ....\n"

else:
	print "\v Recheck the rep1 and rep2 dataframe, aren't of equal sizes....\n"

tf_list1 = rep1_df["Experiment_target"].tolist()
tf_list2 = rep2_df["Experiment_target"].tolist()

### Check if multiple replicates i.e more than 2 present, if yes, then append "_2" as suffix:
new_tf_list1 = []
for each in tf_list1:
	if new_tf_list1.count(each) >=1 :
		print "More than 2 replicates present for TF:", each
		new_tf_list1.append(each+"_2")
 	else:
 		new_tf_list1.append(each)

rep1_df["Experiment_target"] = new_tf_list1

new_tf_list2 = []
for each in tf_list2:
	if new_tf_list2.count(each) >=1 :
		print "More than 2 replicates present for TF:\n", each
		new_tf_list2.append(each+"_2")
 	else:
 		new_tf_list2.append(each)

rep2_df["Experiment_target"] = new_tf_list2


### Create the dummy ids for the File_accession, and Controlled_by id with SLs#:
reps_accession_list = final_read_df["File_accession"].tolist() 
uniq_reps_accession_list = set(reps_accession_list) # unique element list with set()

control_list = []
for each_id in final_read_df["Controlled_by"].tolist():
	if each_id not in uniq_reps_accession_list:
		control_list.append(each_id)

uniq_control_accesion_list = set(control_list)
total_accesion_list = list(uniq_reps_accession_list.union(uniq_control_accesion_list))


### Make the dictionary of accession key and dummy SL values for flexible retrieval:
#SL_list = ["SL000" + str(i) for i in range(1,len(total_accesion_list)+1)] #[map(int, x) for x in T1]
SL_list = ["SL" + str(i) for i in range(SL_start_num, SL_start_num+len(total_accesion_list))]
SL_accession_dict = dict(zip(total_accesion_list, SL_list)) 

"""Using unique total_accession_list with SL# to make SL_directory, and then create a link of fastq files"""
for key,value in SL_accession_dict.iteritems():
	fastq_id = key
	fastq_file = fastq_id + ".fastq.gz"
	if fastq_file not in downloaded_file_list:
		print "\nDownloading the needed fastq....... ::\n", each
		os.environ["fastq_id"] = fastq_id
		os.environ["output_dir"] = output_dir
		os.system("wget --directory-prefix=${output_dir} https://www.encodeproject.org/files/${fastq_id}/@@download/${fastq_id}.fastq.gz")
	
	"""Making the SLs dir and creating the links"""
	#print fastq_file, key
	os.environ["fastq_file"] = fastq_file
	os.environ["SL_dir"] = value
	os.environ["output_dir"] = output_dir
	os.system("if [[ ! -d $output_dir/$SL_dir ]];then mkdir $output_dir/$SL_dir; fi")
	os.system("ln -fs $output_dir/$fastq_file $output_dir/$SL_dir")


### Indexing, or listing the SL values to each accession id:
"""For rep1"""
rep1_SL_list = []
rep1_SL_accesion_list = []

for each in rep1_df["File_accession"].tolist():
	rep1_SL_list.append(SL_accession_dict[each])
	rep1_SL_accesion_list.append((each, SL_accession_dict[each]))

control1_SL_list = []
control1_SL_accesion_list = []
for each in rep1_df["Controlled_by"].tolist():
	control1_SL_list.append(SL_accession_dict[each])
	control1_SL_accesion_list.append((each, SL_accession_dict[each]))

"""For rep2"""
rep2_SL_list = []
rep2_SL_accesion_list = []

for each in rep2_df["File_accession"].tolist():
	rep2_SL_list.append(SL_accession_dict[each])
	rep2_SL_accesion_list.append((each, SL_accession_dict[each]))

control2_SL_list = []
control2_SL_accesion_list = []
for each in rep2_df["Controlled_by"].tolist():
	control2_SL_list.append(SL_accession_dict[each])
	control2_SL_accesion_list.append((each, SL_accession_dict[each]))

### Assigning the SL values, or list to dataframes:
"""For rep1"""
rep1_df["rep_id"] = rep1_SL_list
rep1_df["control_id"] = control1_SL_list
#select_cols = [u'File_accession', u'Controlled_by', u'rep_id', u'control_id', u'Experiment_target',  u'Experiment accession', u'Lab']
select_cols = [u'File_accession', u'Controlled_by', u'rep_id', u'control_id', u'Experiment_target']
rep1_df = rep1_df.loc[:,select_cols]
rep1_df.columns = map(lambda x : x + "_1", rep1_df.columns)

"""For rep2"""
rep2_df["rep_id"] = rep2_SL_list
rep2_df["control_id"] = control2_SL_list
select_cols = [u'File_accession', u'Controlled_by', u'rep_id', u'control_id', u'Experiment_target']
#select_cols = [u'File_accession', u'Controlled_by', u'rep_id', u'control_id', u'Experiment_target', u'Experiment accession', u'Lab']
rep2_df = rep2_df.loc[:,select_cols]
rep2_df.columns = map(lambda x : x + "_2", rep2_df.columns)


final_rep_df = pd.concat([rep1_df.reset_index(), rep2_df.reset_index()], axis=1)
final_rep_df = final_rep_df.drop(["index"], axis=1)
final_rep_df.to_csv(join(output_dir,"final_fastq_download_metadata.txt"), sep="\t", header=True, index=True)

### Save sub-files for ENCODE pipeline processing:
final_rep_df["rep_id_1"].to_csv(join(output_dir, batch_name+"_tf_rep1.txt"), sep="\t", header=None, index=False)
final_rep_df["rep_id_2"].to_csv(join(output_dir, batch_name+"_tf_rep2.txt"), sep="\t", header=None, index=False)
final_rep_df["control_id_1"].to_csv(join(output_dir, batch_name+"_tf_control1.txt"), sep="\t", header=None, index=False)
final_rep_df["control_id_2"].to_csv(join(output_dir, batch_name+"_tf_control2.txt"), sep="\t", header=None, index=False)


