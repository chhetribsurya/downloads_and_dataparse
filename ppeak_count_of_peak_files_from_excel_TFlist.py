#!/usr/bin/env python

import glob
import os
import subprocess
import re
from os.path import join
from os.path import basename

#bsub_options = "bsub -We 10:00 -q normal -n 2 -R span[hosts=1]"

dir_path = os.path.expanduser("~/for_chris/batch_I/idr_passed_peaks_total/SL*")
file_list = glob.glob(dir_path)

#output_dir = os.path.expanduser("~/for_chris/batch_I/chromHMM_overlap_piechart_total")
#if not os.path.exists(output_dir):
#    os.makedirs(output_dir)

excel_file = os.path.expanduser("~/for_chris/batch_I/Encode_full_hepg2_datasets.xlsx")
read_file = pd.read_excel(excel_file, sep="\t", header=True)
print read_file.head()

for
        for each_tf_file in file_list:
            tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_tf_file))[0]
            # general environment variables, preparing for pipeline run:   
            os.environ["peak_file_full_path"] = each_tf_file
            #os.environ["output_dir"] = output_dir 
            #os.environ["piechart_calling_script"] = "/gpfs/gpfs1/home/schhetri/for_encode/spp/analysis_scripts/chromHMM_piechart_hepg2_total.py" 
            output_peakcount = subprocess.check_output('wc -l $peak_file_full_path | cut -f1 -d " " ', shell=True)
            #num_peaks = os.system('wc -l $peak_file_full_path | cut -f1 -d " "')
            print tf_name, ":", output_peakcount