import re
import os
import pandas as pd
import pybedtools
import glob
import numpy
from os.path import join
from os.path import splitext
from os.path import basename

input_file = "download_snyder_fastq_metadata.tsv"
output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_fastqs/snyder"
lab_name = "Snyder" #For all lab; lab_name = ""

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

read_file = pd.read_csv(input_file, sep="\t")
imp_cols = ["File accession", "Experiment accession", "Controlled by", "Biological replicate(s)", "Experiment target", "Biosample term name", \
			"Read length", "Lab", "File download URL",  "File Status",  "Audit NOT_COMPLIANT", "Audit INTERNAL_ACTION", 'Audit WARNING', "md5sum"]
read_file = read_file.loc[:, imp_cols]
select_cols = ["File accession", "Controlled by", "Experiment target", "Biosample term name", "Biological replicate(s)", "Read length", "Experiment accession", "Lab"]
read_df = read_file.loc[:, select_cols]
read_df.columns = read_df.columns.str.replace(" ", "_")

#control_list = read_df['Controlled by'].dropna().replace(r".*files\/(.*)\/", r"\1", regex=True).unique().tolist()
control_list = read_df['Controlled_by'].dropna().replace(r".*files\/([0-9A-Za-z]+)\/*", r"\1", regex=True).unique().tolist()
reps_list = read_df['File_accession'].unique().tolist()

### Download the controls if not available with in the fastq:
for each in control_list:
	if each not in reps_list:
		print "Control Fastq missing:", each+".fastq.gz"
		downloaded_file_list = [ basename(each_fastq) for each_fastq in glob.glob(join(output_dir,"*fastq.gz")) ]
		#check if the fastq was already downloaded:
		if each+".fastq.gz" not in downloaded_file_list:
			print "Downloading control fastq.....:", each
			os.environ["fastq_id"] = each+".fastq.gz"
			os.environ["output_dir"] = "/gpfs1/gpfs1/"
			os.system("wget --directory-prefix=$output_dir https://www.encodeproject.org/files/$fastq_id/@@download/$fastq_id.fastq.gz")
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

pd.concat([rep1_df.reset_index(),rep2_df.reset_index()], axis=1)


### Even or odd slicing of the dataframes: L = [1,2,3,4,5,6,7,8,9,10]; L[start:stop:step]; # prints 1 3 5 7 9
# rep_1_df = final_read_df.iloc[::2] #slicing starting from 0 to all the way, with step 2 i.e even slicing 0,2,4,6...
# rep_2_df = final_read_df.iloc[1::2] #slicing starting from 1 to all the way, with step 2 i.e odd slicing 1,3,5,7...

# mylist = [1,2,3,4,5,6,7,8,9,10]
# for i in mylist[::2]:
#     print i,
# # prints 1 3 5 7 9

# for i in mylist[1::2]:
#     print i,
# # prints 2 4 6 8 10

# Where the first digit is the starting index (defaults to beginning of list or 0), 
# 2nd is ending slice index (defaults to end of list), and the third digit is the offset or step.



# rep1_df["rep1_id"] = [ "SL0000" + str(i) for i in range(1,rep1_df.shape[0]*2,2) ]
# rep1_df["control1_id"] = [ "SL0000" + str(i) for i in range((rep1_df.shape[0]*2 + 1), rep1_df.shape[0]*4, 2) ]

# rep2_df["rep2_id"] = [ "SL0000" + str(i) for i in range(2,rep1_df.shape[0]*2,2) ]
# rep2_df["control2_id"] = [ "SL0000" + str(i) for i in range((rep1_df.shape[0]*2 + 2), rep1_df.shape[0]*4, 2) ]

# ## Even or odd slicing of the dataframes: L = [1,2,3,4,5,6,7,8,9,10]; L[start:stop:step]; # prints 1 3 5 7 9
# rep_1_df = final_read_df.iloc[::2] #slicing starting from 0 to all the way, with step 2 i.e even slicing 0,2,4,6...
# rep_2_df = final_read_df.iloc[1::2] #slicing starting from 1 to all the way, with step 2 i.e odd slicing 1,3,5,7...




#############################

#### Alternative approach to tackle the problem of extracting accurate biological replicates, and edited for more than 2
#### bilogical replicates:

#############################


input_file = "download_snyder_fastq_metadata.tsv"
output_dir = "/gpfs/gpfs1/home/schhetri/for_SLs/batch_III_libraries/download_fastqs/snyder"
lab_name = "Snyder" #For all lab; lab_name = ""

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

read_file = pd.read_csv(input_file, sep="\t")
imp_cols = ["File accession", "Experiment accession", "Controlled by", "Biological replicate(s)", "Experiment target", "Biosample term name", \
			"Read length", "Lab", "File download URL",  "File Status",  "Audit NOT_COMPLIANT", "Audit INTERNAL_ACTION", 'Audit WARNING', "md5sum"]
read_file = read_file.loc[:, imp_cols]
select_cols = ["File accession", "Controlled by", "Experiment target", "Biosample term name", "Biological replicate(s)", "Read length", "Experiment accession", "Lab"]
read_df = read_file.loc[:, select_cols]
read_df.columns = read_df.columns.str.replace(" ", "_")

#control_list = read_df['Controlled by'].dropna().replace(r".*files\/(.*)\/", r"\1", regex=True).unique().tolist()
control_list = read_df['Controlled_by'].dropna().replace(r".*files\/([0-9A-Za-z]+)\/*", r"\1", regex=True).unique().tolist()
reps_list = read_df['File_accession'].unique().tolist()

### Download the controls if not available with in the fastq:
for each in control_list:
	if each not in reps_list:
		print "Control Fastq missing:", each+".fastq.gz"
		downloaded_file_list = [ basename(each_fastq) for each_fastq in glob.glob(join(output_dir,"*fastq.gz")) ]
		#check if the fastq was already downloaded:
		if each+".fastq.gz" not in downloaded_file_list:
			print "Downloading control fastq.....:", each
			os.environ["fastq_id"] = each+".fastq.gz"
			os.environ["output_dir"] = "/gpfs1/gpfs1/"
			os.system("wget --directory-prefix=$output_dir https://www.encodeproject.org/files/$fastq_id/@@download/$fastq_id.fastq.gz")
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

final_read_df = read_df.iloc[rep_idx]

### Select the samples with replicates:
rep1 = final_read_df[final_read_df["Biological_replicate(s)"] == 1]
rep2 = final_read_df[final_read_df["Biological_replicate(s)"] == 2]
tf_list1 = rep1["Experiment_target"].tolist()
tf_list2 = rep2["Experiment_target"].tolist()

new_tf_list1 = []
for each in tf_list1:
	if new_tf_list1.count(each) >=1 :
		new_tf_list1.append(each+"_2")
 	else:
 		new_tf_list1.append(each)

rep1["Experiment_target"] = new_tf_list1


new_tf_list2 = []
for each in tf_list2:
	if new_tf_list2.count(each) >=1 :
		new_tf_list2.append(each+"_2")
 	else:
 		new_tf_list2.append(each)

rep2["Experiment_target"] = new_tf_list2


#############################








In [248]: rep2.head()
Out[248]:
   File_accession Controlled_by Experiment_target Biosample_term_name  \
0     ENCFF002EAS   ENCFF002ECQ      NFE2L2-human               HepG2
3     ENCFF002EDS   ENCFF002ECQ       HCFC1-human               HepG2
5     ENCFF000XTZ   ENCFF000XTF        MAFF-human               HepG2
10    ENCFF000XQN   ENCFF000XTF       CEBPB-human               HepG2
12    ENCFF000XWR   ENCFF000XTF      POLR2A-human               HepG2

    Biological_replicate(s)  Read_length Experiment_accession  \
0                         2           36          ENCSR488EES
3                         2           36          ENCSR529JYA
5                         2           36          ENCSR000EEC
10                        2           36          ENCSR000EEE
12                        2           36          ENCSR000EEM

                         Lab
0   Michael Snyder, Stanford
3   Michael Snyder, Stanford
5   Michael Snyder, Stanford
10  Michael Snyder, Stanford
12  Michael Snyder, Stanford

In [249]: rep1.head()
Out[249]:
   File_accession Controlled_by Experiment_target Biosample_term_name  \
1     ENCFF002EAY   ENCFF002ECU      NFE2L2-human               HepG2
2     ENCFF002EDN   ENCFF002ECU       HCFC1-human               HepG2
6     ENCFF000XUA   ENCFF000XTF        MAFF-human               HepG2
9     ENCFF000XQM   ENCFF000XTF       CEBPB-human               HepG2
11    ENCFF000XWJ   ENCFF000XTF      POLR2A-human               HepG2

    Biological_replicate(s)  Read_length Experiment_accession  \
1                         1           36          ENCSR488EES
2                         1           36          ENCSR529JYA
6                         1           36          ENCSR000EEC
9                         1           36          ENCSR000EEE
11                        1           36          ENCSR000EEM

                         Lab
1   Michael Snyder, Stanford
2   Michael Snyder, Stanford
6   Michael Snyder, Stanford
9   Michael Snyder, Stanford
11  Michael Snyder, Stanford

In [250]: rep1.tail()
Out[250]:
    File_accession Controlled_by Experiment_target Biosample_term_name  \
93     ENCFF000XUF   ENCFF000XTF         MAX-human               HepG2
95     ENCFF000XQP   ENCFF000XTF        CHD2-human               HepG2
97     ENCFF000XPL   ENCFF000XTF       BRCA1-human               HepG2
99     ENCFF000XXY   ENCFF000XTF        SMC3-human               HepG2
102    ENCFF000XRO   ENCFF000XQR       ESRRA-human               HepG2

     Biological_replicate(s)  Read_length Experiment_accession  \
93                         1           36          ENCSR000EDS
95                         1           36          ENCSR000EED
97                         1           36          ENCSR000EDY
99                         1           36          ENCSR000EDW
102                        1           27          ENCSR000EEW

                          Lab
93   Michael Snyder, Stanford
95   Michael Snyder, Stanford
97   Michael Snyder, Stanford
99   Michael Snyder, Stanford
102  Michael Snyder, Stanford

In [251]: rep2.tail()
Out[251]:
    File_accession Controlled_by Experiment_target Biosample_term_name  \
94     ENCFF000XUH   ENCFF000XTF         MAX-human               HepG2
96     ENCFF000XQQ   ENCFF000XTF        CHD2-human               HepG2
98     ENCFF000XPS   ENCFF000XTF       BRCA1-human               HepG2
100    ENCFF000XYC   ENCFF000XTF        SMC3-human               HepG2
101    ENCFF000XRN   ENCFF000XQS       ESRRA-human               HepG2

     Biological_replicate(s)  Read_length Experiment_accession  \
94                         2           36          ENCSR000EDS
96                         2           36          ENCSR000EED
98                         2           36          ENCSR000EDY
100                        2           36          ENCSR000EDW
101                        2           27          ENCSR000EEW

                          Lab
94   Michael Snyder, Stanford
96   Michael Snyder, Stanford
98   Michael Snyder, Stanford
100  Michael Snyder, Stanford
101  Michael Snyder, Stanford

In [252]: rep2.shape
Out[252]: (45, 8)

In [253]: rep1.shape
Out[253]: (45, 8)



In [258]: new_tf_list = []

In [259]: for each in tf_list:
   .....:     if new_tf_list.count(each) >=1 :
   .....:         new_tf_list.append(each+"_2")
   .....:     else:
   .....:         new_tf_list.append(each)
   .....:

In [260]: new_tf_list
Out[260]:
['NFE2L2-human',
 'HCFC1-human',
 'MAFF-human',
 'CEBPB-human',
 'POLR2A-human',
 'TBP-human',
 'NR3C1-human',
 'RAD21-human',
 'HNF4A-human',
 'POLR2AphosphoS2-human',
 'USF2-human',
 'MAFK-human',
 'IRF3-human',
 'CEBPZ-human',
 'SUZ12-human',
 'MXI1-human',
 'EP300-human',
 'RCOR1-human',
 'KAT2B-human',
 'GTF2F1-human',
 'RFX5-human',
 'JUN-human',
 'FOS-human',
 'ZNF384-human',
 'TBL1XR1-human',
 'ZNF143-human',
 'MAFK-human_2',
 'POLR2A-human_2',
 'MAZ-human',
 'JUND-human',
 'CEBPB-human_2',
 'CUX1-human',
 'BACH1-human',
 'SIN3B-human',
 'ARID3A-human',
 'BHLHE40-human',
 'PPARGC1A-human',
 'HSF1-human',
 'NRF1-human',
 'SREBF1-human',
 'MAX-human',
 'CHD2-human',
 'BRCA1-human',
 'SMC3-human',
 'ESRRA-human']





### Edit all the tfs with more than 2 replicate counts with additional suffix:
new_tf_list = []
for item in tf_list:
	if new_tf_list.count(item) >= 2:
		new_tf_list.append(item + "_2")
	else:
		new_tf_list.append(item)

final_rep_df["Experiment_target"] = new_tf_list

In [181]: new_tf_list
Out[181]:
['NFE2L2-human',
 'NFE2L2-human',
 'HCFC1-human',
 'HCFC1-human',
 'MAFF-human',
 'MAFF-human',
 'CEBPB-human',
 'CEBPB-human',
 'POLR2A-human',
 'POLR2A-human',
 'TBP-human',
 'TBP-human',
 'NR3C1-human',
 'NR3C1-human',
 'RAD21-human',
 'RAD21-human',
 'HNF4A-human',
 'HNF4A-human',
 'POLR2AphosphoS2-human',
 'POLR2AphosphoS2-human',
 'USF2-human',
 'USF2-human',
 'MAFK-human',
 'MAFK-human',
 'IRF3-human',
 'IRF3-human',
 'CEBPZ-human',
 'CEBPZ-human',
 'SUZ12-human',
 'SUZ12-human',
 'MXI1-human',
 'MXI1-human',
 'MXI1-human_2',
 'EP300-human',
 'EP300-human',
 'RCOR1-human',
 'RCOR1-human',
 'KAT2B-human',
 'KAT2B-human',
 'GTF2F1-human',
 'GTF2F1-human',
 'ZKSCAN1-human',
 'ZKSCAN1-human',
 'RFX5-human',
 'RFX5-human',
 'JUN-human',
 'JUN-human',
 'FOS-human',
 'FOS-human',
 'ZNF384-human',
 'ZNF384-human',
 'TBL1XR1-human',
 'TBL1XR1-human',
 'ZNF143-human',
 'ZNF143-human',
 'MAFK-human_2',
 'MAFK-human_2',
 'POLR2A-human_2',
 'POLR2A-human_2',
 'MAZ-human',
 'MAZ-human',
 'JUND-human',
 'JUND-human',
 'CEBPB-human_2',
 'CEBPB-human_2',
 'CUX1-human',
 'CUX1-human',
 'BACH1-human',
 'BACH1-human',
 'SIN3B-human',
 'SIN3B-human',
 'ARID3A-human',
 'ARID3A-human',
 'BHLHE40-human',
 'BHLHE40-human',
 'PPARGC1A-human',
 'PPARGC1A-human',
 'HSF1-human',
 'HSF1-human',
 'NRF1-human',
 'NRF1-human',
 'SREBF1-human',
 'SREBF1-human',
 'MAX-human',
 'MAX-human',
 'CHD2-human',
 'CHD2-human',
 'BRCA1-human',
 'BRCA1-human',
 'SMC3-human',
 'SMC3-human',
 'ESRRA-human',
 'ESRRA-human']




tuple_list = []
for i,item in enumerate(final_tf_list):
	if final_tf_list.count(item) > 2:
		print item, i
		tuple_list.append((item, i))







more_then_2rep_idx = []
for i,item in enumerate(tf_name_list):
	if tf_name_list.count(item) > 2:
		print(item)
		more_then_2rep_idx.append(i)

more_rep_df = read_df.iloc[more_then_2rep_idx].loc[:,["File accession",'Controlled by',"Experiment target"]]

final_df = pd.concat([2rep_df,more_rep_df], ignore_index=True)







sorted_file = read_df.sort_values(["Experiment target"])
tf_name_list = sorted_file.loc[:,["Experiment target"]]

single_rep_idx = []
for i,item in enumerate(tf_name_list):
	if tf_name_list.count(item) < 2:
		print(item)
		single_rep_idx.append(i)

single_rep_df = sorted_file.iloc[single_rep_idx].loc[:,["File accession",'Controlled by',"Experiment target"]]


2rep_idx = []
for i,item in enumerate(tf_name_list):
	if tf_name_list.count(item) == 2:
		print(item)
		2rep_idx.append(i)

2rep_df = sorted_file.iloc[2rep_idx].loc[:,["File accession",'Controlled by',"Experiment target"]]



more_then_2rep_idx = []
for i,item in enumerate(tf_name_list):
	if tf_name_list.count(item) > 2:
		print(item)
		more_then_2rep_idx.append(i)

more_rep_df = sorted_file.iloc[more_then_2rep_idx].loc[:,["File accession",'Controlled by',"Experiment target"]]


for idx,tf in enumerate(rep_list):
	each_tf_list = []
	if rep_list.count(tf) == 4:
		each_tf_list.append(idx)
		print each, i



CEBPB-human 6
CEBPB-human 7
POLR2A-human 8
POLR2A-human 9
MAFK-human 22
MAFK-human 23
MAFK-human 52
MAFK-human 53
POLR2A-human 54
POLR2A-human 55
CEBPB-human 60
CEBPB-human 61


In [38]: sorted_file.iloc[more_then_2rep_idx].loc[:,["File accession",'Controlled by',"Experiment target"]]
Out[38]:
   File accession        Controlled by Experiment target
9     ENCFF000XQM  /files/ENCFF000XTF/       CEBPB-human
10    ENCFF000XQN  /files/ENCFF000XTF/       CEBPB-human
11    ENCFF000XWJ  /files/ENCFF000XTF/      POLR2A-human
12    ENCFF000XWR  /files/ENCFF000XTF/      POLR2A-human
25    ENCFF000XUK  /files/ENCFF000XTF/        MAFK-human
26    ENCFF000XUL  /files/ENCFF000XTF/        MAFK-human
33    ENCFF000XUX  /files/ENCFF000XTE/        MXI1-human
34    ENCFF000XVA  /files/ENCFF000XTE/        MXI1-human
35    ENCFF000XVC  /files/ENCFF000XTE/        MXI1-human
36    ENCFF000XRH                  NaN     Control-human
37    ENCFF000XRI                  NaN     Control-human
42    ENCFF000XQR                  NaN     Control-human
43    ENCFF000XQS                  NaN     Control-human
60    ENCFF002ECQ                  NaN     Control-human
61    ENCFF002ECU                  NaN     Control-human
64    ENCFF000XUU  /files/ENCFF000XTF/        MAFK-human
65    ENCFF000XVK  /files/ENCFF000XTF/        MAFK-human
66    ENCFF000XWF  /files/ENCFF000XQS/      POLR2A-human
67    ENCFF000XWI  /files/ENCFF000XQR/      POLR2A-human
72    ENCFF000XPR  /files/ENCFF000XQR/       CEBPB-human
73    ENCFF000XPT  /files/ENCFF000XQS/       CEBPB-human
90    ENCFF000XYA                  NaN      SREBF1-human
91    ENCFF000XYB  /files/ENCFF000XRH/      SREBF1-human
92    ENCFF000XYE  /files/ENCFF000XRI/      SREBF1-human





In [39]: multiple_reps_df = read_file.iloc[idx_list].loc[:,["File accession",'Controlled by',"Experiment target"]]


In [41]: multiple_reps_df.groupby('Experiment target').nth(1)
Out[41]:
                         Controlled by File accession
Experiment target
CEBPB-human        /files/ENCFF000XTF/    ENCFF000XQN
Control-human                      NaN    ENCFF000XRI
MAFK-human         /files/ENCFF000XTF/    ENCFF000XUL
MXI1-human         /files/ENCFF000XTE/    ENCFF000XVA
POLR2A-human       /files/ENCFF000XTF/    ENCFF000XWR
SREBF1-human       /files/ENCFF000XRH/    ENCFF000XYB


### Limit the no. of occurences in pandas:( Pandas dataframe get first row of each group )
In [42]: multiple_reps_df.groupby('Experiment target').nth([1,2])
Out[42]:
                         Controlled by File accession
Experiment target
CEBPB-human        /files/ENCFF000XTF/    ENCFF000XQN
CEBPB-human        /files/ENCFF000XQR/    ENCFF000XPR
Control-human                      NaN    ENCFF000XRI
Control-human                      NaN    ENCFF000XQR
MAFK-human         /files/ENCFF000XTF/    ENCFF000XUL
MAFK-human         /files/ENCFF000XTF/    ENCFF000XUU
MXI1-human         /files/ENCFF000XTE/    ENCFF000XVA
MXI1-human         /files/ENCFF000XTE/    ENCFF000XVC
POLR2A-human       /files/ENCFF000XTF/    ENCFF000XWR
POLR2A-human       /files/ENCFF000XQS/    ENCFF000XWF
SREBF1-human       /files/ENCFF000XRH/    ENCFF000XYB
SREBF1-human       /files/ENCFF000XRI/    ENCFF000XYE








In [41]: read_file[read_file['Controlled by'].isnull()]["Experiment target"]
Out[41]:
4     rabbit-IgG-control-human
7       goat-IgG-control-human
8       goat-IgG-control-human
36               Control-human
37               Control-human
42               Control-human
43               Control-human
60               Control-human
61               Control-human
90                SREBF1-human
Name: Experiment target, dtype: object

In [42]: control_df = read_file[read_file['Controlled by'].isnull()]["Experiment target"]

In [43]: read_file[~read_file['Controlled by'].isnull()].loc[:,["File accession",'Controlled by',"Experiment target"]]
Out[43]:
    File accession        Controlled by      Experiment target
0      ENCFF002EAS  /files/ENCFF002ECQ/           NFE2L2-human
1      ENCFF002EAY  /files/ENCFF002ECU/           NFE2L2-human
2      ENCFF002EDN  /files/ENCFF002ECU/            HCFC1-human
3      ENCFF002EDS  /files/ENCFF002ECQ/            HCFC1-human
5      ENCFF000XTZ  /files/ENCFF000XTF/             MAFF-human
6      ENCFF000XUA  /files/ENCFF000XTF/             MAFF-human
9      ENCFF000XQM  /files/ENCFF000XTF/            CEBPB-human
10     ENCFF000XQN  /files/ENCFF000XTF/            CEBPB-human
11     ENCFF000XWJ  /files/ENCFF000XTF/           POLR2A-human
12     ENCFF000XWR  /files/ENCFF000XTF/           POLR2A-human
13     ENCFF000XZI  /files/ENCFF000XTF/              TBP-human
14     ENCFF000XZL  /files/ENCFF000XTF/              TBP-human
15     ENCFF000XRU  /files/ENCFF000XQS/            NR3C1-human
16     ENCFF000XRV  /files/ENCFF000XQR/            NR3C1-human
17     ENCFF000XXK  /files/ENCFF000XTF/            RAD21-human
18     ENCFF000XXL  /files/ENCFF000XTF/            RAD21-human
19     ENCFF000XRT  /files/ENCFF000XQR/            HNF4A-human
20     ENCFF000XSM  /files/ENCFF000XQS/            HNF4A-human
21     ENCFF000XXB  /files/ENCFF000XTF/  POLR2AphosphoS2-human
22     ENCFF000XXC  /files/ENCFF000XTF/  POLR2AphosphoS2-human
23     ENCFF000XZT  /files/ENCFF000XTF/             USF2-human
24     ENCFF000XZV  /files/ENCFF000XTF/             USF2-human
25     ENCFF000XUK  /files/ENCFF000XTF/             MAFK-human
26     ENCFF000XUL  /files/ENCFF000XTF/             MAFK-human
27     ENCFF000XTN  /files/ENCFF000XTF/             IRF3-human
28     ENCFF000XTO  /files/ENCFF000XTF/             IRF3-human
29     ENCFF000XPZ  /files/ENCFF000XTF/            CEBPZ-human
30     ENCFF000XQB  /files/ENCFF000XTF/            CEBPZ-human
31     ENCFF002EEA  /files/ENCFF002ECU/            SUZ12-human
32     ENCFF002EEF  /files/ENCFF002ECQ/            SUZ12-human
..             ...                  ...                    ...
72     ENCFF000XPR  /files/ENCFF000XQR/            CEBPB-human
73     ENCFF000XPT  /files/ENCFF000XQS/            CEBPB-human
74     ENCFF002EDU  /files/ENCFF002ECU/             CUX1-human
75     ENCFF002EDX  /files/ENCFF002ECQ/             CUX1-human
76     ENCFF002EDO  /files/ENCFF002ECQ/            BACH1-human
77     ENCFF002EDR  /files/ENCFF002ECU/            BACH1-human
78     ENCFF002EDW  /files/ENCFF002ECU/            SIN3B-human
79     ENCFF002EDY  /files/ENCFF002ECQ/            SIN3B-human
80     ENCFF000XOS  /files/ENCFF000XTF/           ARID3A-human
81     ENCFF000XOU  /files/ENCFF000XTF/           ARID3A-human
82     ENCFF000XPA  /files/ENCFF000XTF/          BHLHE40-human
83     ENCFF000XPC  /files/ENCFF000XTF/          BHLHE40-human
84     ENCFF000XVQ  /files/ENCFF000XQR/         PPARGC1A-human
85     ENCFF000XWG  /files/ENCFF000XQS/         PPARGC1A-human
86     ENCFF000XSK  /files/ENCFF000XQS/             HSF1-human
87     ENCFF000XSL  /files/ENCFF000XQR/             HSF1-human
88     ENCFF000XVP  /files/ENCFF000XTF/             NRF1-human
89     ENCFF000XVR  /files/ENCFF000XTF/             NRF1-human
91     ENCFF000XYB  /files/ENCFF000XRH/           SREBF1-human
92     ENCFF000XYE  /files/ENCFF000XRI/           SREBF1-human
93     ENCFF000XUF  /files/ENCFF000XTF/              MAX-human
94     ENCFF000XUH  /files/ENCFF000XTF/              MAX-human
95     ENCFF000XQP  /files/ENCFF000XTF/             CHD2-human
96     ENCFF000XQQ  /files/ENCFF000XTF/             CHD2-human
97     ENCFF000XPL  /files/ENCFF000XTF/            BRCA1-human
98     ENCFF000XPS  /files/ENCFF000XTF/            BRCA1-human
99     ENCFF000XXY  /files/ENCFF000XTF/             SMC3-human
100    ENCFF000XYC  /files/ENCFF000XTF/             SMC3-human
101    ENCFF000XRN  /files/ENCFF000XQS/            ESRRA-human
102    ENCFF000XRO  /files/ENCFF000XQR/            ESRRA-human




In [54]: read_file['Controlled by'].str.replace("/files/","").replace("/","").unique().tolist()
Out[54]:
['ENCFF002ECQ/',
 'ENCFF002ECU/',
 nan,
 'ENCFF000XTF/',
 'ENCFF000XQS/',
 'ENCFF000XQR/',
 'ENCFF000XTE/',
 'ENCFF000XRH/',
 'ENCFF000XRI/']

In [55]: read_file['Controlled by'].str.replace("/files/","").str.replace("/","").unique().tolist()
Out[55]:
['ENCFF002ECQ',
 'ENCFF002ECU',
 nan,
 'ENCFF000XTF',
 'ENCFF000XQS',
 'ENCFF000XQR',
 'ENCFF000XTE',
 'ENCFF000XRH',
 'ENCFF000XRI']

In [56]: control_list = read_file['Controlled by'].dropna().str.replace("/files/","").str.replace("/","").unique().tolist()


In [12]: reps_list=read_file['File accession'].unique().tolist()

In [57]: read_file[read_file['File accession'].isin(control_list)]


In [63]: for each_ctrl in control_list:
    ...:     if each_ctrl not in reps_list:
    ...:         print(each_ctrl)
    ...:
nan
ENCFF000XTE


In [68]: for each_ctrl in control_list:
    ...:     if each_ctrl in reps_list:
    ...:         print(each_ctrl)
    ...:
ENCFF002ECQ
ENCFF002ECU
ENCFF000XTF
ENCFF000XQS
ENCFF000XQR
ENCFF000XRH
ENCFF000XRI





##############################
##############################
##############################


In [73]: read_file['Controlled by'].dropna().str.replace("/files/","").str.replace("/","").unique().tolist()
Out[73]:
['ENCFF002ECQ',
 'ENCFF002ECU',
 'ENCFF000XTF',
 'ENCFF000XQS',
 'ENCFF000XQR',
 'ENCFF000XTE',
 'ENCFF000XRH',
 'ENCFF000XRI']


 for each in control_list:
	if each not in reps_list:
		print(each)
		os.environ["fastq_id"]=each
		os.environ["output_dir"]="/gpfs1/gpfs1/"
		os.system("wget --directory-prefix=$output_dir https://www.encodeproject.org/files/$fastq_id/@@download/$fastq_id.fastq.gz")


In [113]: total_accesion_list
Out[113]:
{'ENCFF000XOS',
 'ENCFF000XOU',
 'ENCFF000XPA',
 'ENCFF000XPC',
 'ENCFF000XPL',
 'ENCFF000XPR',
 'ENCFF000XPS',
 'ENCFF000XPT',
 'ENCFF000XPZ',
 'ENCFF000XQB',
 'ENCFF000XQM',
 'ENCFF000XQN',
 'ENCFF000XQP',
 'ENCFF000XQQ',
 'ENCFF000XQR',
 'ENCFF000XQS',
 'ENCFF000XQT',
 'ENCFF000XQU',
 'ENCFF000XQW',
 'ENCFF000XQY',
 'ENCFF000XRH',
 'ENCFF000XRI',
 'ENCFF000XRN',
 'ENCFF000XRO',
 'ENCFF000XRT',
 'ENCFF000XRU',
 'ENCFF000XRV',
 'ENCFF000XSK',
 'ENCFF000XSL',
 'ENCFF000XSM',
 'ENCFF000XTE',
 'ENCFF000XTF',
 'ENCFF000XTN',
 'ENCFF000XTO',
 'ENCFF000XTQ',
 'ENCFF000XTR',
 'ENCFF000XTZ',
 'ENCFF000XUA',
 'ENCFF000XUF',
 'ENCFF000XUH',
 'ENCFF000XUK',
 'ENCFF000XUL',
 'ENCFF000XUN',
 'ENCFF000XUP',
 'ENCFF000XUU',
 'ENCFF000XUX',
 'ENCFF000XVA',
 'ENCFF000XVC',
 'ENCFF000XVK',
 'ENCFF000XVP',
 'ENCFF000XVQ',
 'ENCFF000XVR',
 'ENCFF000XVS',
 'ENCFF000XVT',
 'ENCFF000XWF',
 'ENCFF000XWG',
 'ENCFF000XWI',
 'ENCFF000XWJ',
 'ENCFF000XWR',
 'ENCFF000XXB',
 'ENCFF000XXC',
 'ENCFF000XXK',
 'ENCFF000XXL',
 'ENCFF000XXX',
 'ENCFF000XXY',
 'ENCFF000XXZ',
 'ENCFF000XYB',
 'ENCFF000XYC',
 'ENCFF000XYE',
 'ENCFF000XZI',
 'ENCFF000XZL',
 'ENCFF000XZT',
 'ENCFF000XZV',
 'ENCFF001YXD',
 'ENCFF001YXE',
 'ENCFF001YXF',
 'ENCFF001YXH',
 'ENCFF001YXI',
 'ENCFF001YXJ',
 'ENCFF002EAR',
 'ENCFF002EAS',
 'ENCFF002EAV',
 'ENCFF002EAY',
 'ENCFF002ECQ',
 'ENCFF002ECU',
 'ENCFF002EDN',
 'ENCFF002EDO',
 'ENCFF002EDR',
 'ENCFF002EDS',
 'ENCFF002EDT',
 'ENCFF002EDU',
 'ENCFF002EDV',
 'ENCFF002EDW',
 'ENCFF002EDX',
 'ENCFF002EDY',
 'ENCFF002EDZ',
 'ENCFF002EEA',
 'ENCFF002EED',
 'ENCFF002EEF',
 'ENCFF002EEG',
 'ENCFF002EEH'}







