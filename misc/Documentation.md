# Setup Requirements

### Initial setups

In the home directory create the folder *.seqkit* and place [*seqkit.yaml*](https://github.com/ashwini06/seqkit/blob/master/data/seqkit.yaml) (configuration file) inside *.seqkit* folder. 
Edit the root dir path in *seqkit.yaml*

It should point to the path one level up where all the data folders are present

Eg:  `root_dir: "/proj/b2012025/RAW_DATA/ChIP_Data"`


To run seqkit we need some folder re-arrangements for placing the rawdata files. Seqkit considers [this](https://www.evernote.com/shard/s734/res/a0538341-8226-4583-8d3d-559c31a6b476/Seqkit_project_dir.pdf) file structure. 

## Available functions in seqkit
 
`seqkit  --help`

[Running fastqc](#preqc)

[Aligning and running bam-to-bed (default runs bowtie)](#align)

[Peak-calling (either in TF or HM mode: macs2/danpos2)](#Peak-call)

[Post-QC for analysis](#postqc)

[Runing Tophat, Cufflinks](#none)

[Running htSeq-Count] (#none)

<a name="preqc"></a>
### Running Fastqc 

A prior quality check to run a fastqc on all the sample folders present in the project folder. It is important to check the sequency quality, over-represented sequences and duplicate percentages from the fastqc output results.
Fastqc output explaination is summarized [here](http://www4.ncsu.edu/~rosswhet/BIT815/Overview/Week2/FastQC_details.pdf)

`seqkit preqc -p Ascl1_US`


**Inputs**

>| Command | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |

**Outputs**

Creates fastqc folder inside the sample folder. _*.html_ contains the fastqc summarized report.



<a name="align"/></a>
### Aligning and running bam-to-bed (default runs bowtie)

*Picard tools* is incorporated into _seqkit_ to estimate the quality metrics and remove the duplicates from aligned bam files.
However both the bam files _(*_sorted.bam)_ and _(*_rmdup.bam)_ files are present in the _(alignment_*)_ folder.
For further steps in seqkit, duplicates removed bam files are used .

**To run on all samples present in the project folder (-p)** 

`seqkit analysis align -p Ascl1_US --bam_to_bed`

Bed files are input to the peak calling software (macs2/danpos2).
So bed files can be generated while aligning the reads by adding an extra option *--bam_to_bed*.

Also, the above can independently run (in case if we dont to generate *bed files* 
or if we want to run generate *bed files* from already existed *aligned reads*)

**Command-line**

`seqkit analysis align -p Ascl1_US ` (Doesnt generate bed files)

`seqkit analysis bamtobed -p Ascl1_US --slurm` (Generates bed files from already existed files)

**To run on specific samples (s) present in the project folder (-p)**

`seqkit analysis align -p Ascl1_US --bam_to_bed -s Mark_Mash1_s4`

`seqkit analysis align -p Ascl1_US --bam_to_bed -s Mark_input_s3`


**Inputs**

>| Command | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |
| -s/--sample | FILENAME *(optional)*| to run on specific samples inside project folder |
|--bamtobed | *(optional)*|   to generate bed files |


**Outputs files**
<ul>
<li> Creates alignment_(aligner) folder in the project specific folder </li>
<li> Creates seperate folder for _(bam_files)_ and _(bedfiles)_</li>
</ul>

*/project_dir/sample_folder/alignment_bowtie/bam_files* contains:
+ *_bowtie2.log
+ *_sorted.bam.bai
- *_sorted.bam
- *_sorted_rmdup.bam
- *_sorted_rmdup.bam.bai

*/project_dir/sample_folder/alignment_bowtie/bedfiles* contains:
+ *uniq.bed

<a name="Peak-call"/></a>

### Peak-calling (either in TF or HM mode: macs2/danpos2)

For _TF chip_ the macs2 is working well to call the enriched regions.
Whereas for some of the broad marks, macs2 doesn't work well in calling the enriched regions .
However it calls the regions but it is divided into small segmental peaks.
danpos2 might work well so included danpos2 option for HM data.

**Command-line**

`seqkit peakcall -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt`

**Inputs**

>| Command | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |
| -i/--input | FILENAME | Path to tab-delimited text file. It is a 2 column file, where first column should contain sample name (Treatment) and second column should contain sample name (Control) |
|-m/--mode | TF/HM |  TF chip or Histone modification (default is TF) |
|--peak_call| macs2/danpos2|macs2 (applicable for TF and HM) or danpos2(only for HM) (default is HM) |

**Outputs**


<a name="postqc"/></a>
### Post-QC for analysis

DeepTools provides number of quality metrics and provides an estimate for assessing the quality of ChIP.
Incoporated _bamcompare_,_computeMatrix_,_plotHeatmap_ functions from *deepTools*.

**Command-line**

`seqkit postqc -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt --genefile /home/ashwini/scripts/seqkit/data/UCSC_mm10genes_v2.bed`

**Inputs**


