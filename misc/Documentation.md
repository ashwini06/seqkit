# Setup Requirements

### Initial setups

In the home directory create the folder *.seqkit* and place [seqkit.yaml](https://github.com/ashwini06/seqkit/blob/master/data/seqkit.yaml) (configuration file) inside .seqkit folder


Edit the root dir path in seqkit.yaml:
It should point to the path one level up where all the data folders are present
root_dir: "/proj/b2012025/RAW_DATA/ChIP_histone"


To run seqkit we need some folder rearrangements for placing the rawdata files. [Folder structure](https://www.evernote.com/shard/s734/res/a0538341-8226-4583-8d3d-559c31a6b476/Seqkit_project_dir.pdf)
Sample folder are placed under Project folder and inside each sample folder should contain Rawdata folder and sample folder should contain fastq or fastq.gz files.

seqkit  --help
Gives the available methods in seqkit

Aligning and running bam-to-bed (default runs bowtie)
Picard tools is added up to estimate the quality metrics and remove the duplicates from aligned bam files.
However both the bam files *_sorted.bam and *_rmdup.bam files are present in the alignment_(aligner) folder.
For further steps, duplicates removed bam files are used .

To run on all samples present in the project folder (-p)
seqkit analysis align -p Ascl1_US --bam_to_bed

--bam_to_bed (Optional) : Bed files are input to the peak calling software (macs2/danpos2).
So bed files can be generated while aligning the reads by adding extra option --bam_to_bed.

Also, the above can independently run (in case if we dont to generate bed files or if we want to run generate bed files from already existed aligned reads)

seqkit analysis align -p Ascl1_US  (Doesnt generate bed files)
seqkit analysis bamtobed -p Ascl1_US --slurm (Generates bed files from already existed files)

To run on specific samples (s) present in the project folder (-p)
seqkit analysis align -p Ascl1_US --bam_to_bed -s Mark_Mash1_s4
seqkit analysis align -p Ascl1_US --bam_to_bed -s Mark_input_s3


Inputs
-p / --project : project folder name
-s/ --sample : (optional) to run on specific samples inside project folder
--bamtobed : (optional) to generate bed files

Outputs : (by default)
Creates alignment_(aligner) folder in the project specific folder
Creates seperate folder for bam_files and bedfiles

/project_dir/sample_folder/alignment_bowtie/bam_files:
Mark_input_s3_Markinput_bowtie2.log 
Mark_input_s3_Markinput_sorted.bam.bai
Mark_input_s3_Markinput_sorted.bam

/project_dir/sample_folder/alignment_bowtie/bedfiles:



Peak-calling (either in TF or HM mode: macs2/danpos2)
For TF chip the macs2 is working well to call the enriched regions.
Whereas for some of the broad marks, macs2 doesn't work well in calling the enriched regions .
However it calls the regions but it is small segmental peaks.
Other people suggested that danpos2 might work well so included danpos2 option for HM data.

Inputs:

-p / --project : Project folder
-i / --input : Path to tab-delimited text file. It is a 2 column file, where first column should contain sample name (Treatment) and second column should contain sample name (Control)
-m / --mode : TF chip or Histone modification (default is TF)
--peak_call : macs2 (applicable for TF and HM) or danpos2(only for HM) (default is HM)

Outputs:
seqkit peakcall -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt

Post-QC for analysis
DeepTools provides number of quality metrics and provides an estimate for assessing the quality of ChIP.
Incoporated bamcompare,computeMatrix,plotHeatmap functions from deepTools.

Inputs

seqkit postqc -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt --genefile /home/ashwini/scripts/seqkit/data/UCSC_mm10genes_v2.bed 
