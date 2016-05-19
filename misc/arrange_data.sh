#!/bin/bash -l

# Script to rearrange the folder into rawdata folder

# usage : bash arrange_data.sh /proj/b2012025/RAW_DATA/ChIP_histone/lizzy_andersson_ChIP_ESC_wt_57-64sonctrl/rawdata /proj/b2012025/RAW_DATA/ChIP_histone/lizzy_andersson_ChIP_ESC_wt_57-64sonctrl

# dir_structure : proj_folder/sample_folder/Raw_data/*fastq ; 

# $dir should point to the path where all the sample folders are present.

dir=$1
nw_dir=$2

samples=$(ls --color=never $dir);
for sam in $samples;do
nm=${dir}/${sam}
raw_dir=${dir}/${sam}'/Rawdata'
if [ ! -d  ${raw_dir} ]; then mkdir -p ${raw_dir};fi
mv -t $raw_dir $nm/*.*
mv $nm $nw_dir
done
 






