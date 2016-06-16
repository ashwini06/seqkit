# Setup Requirements

### Initial setup

In the home directory create the folder *.seqkit* and place [*seqkit.yaml*](https://github.com/ashwini06/seqkit/blob/master/data/seqkit.yaml) (configuration file) inside *.seqkit* folder. 

`mkdir ~/.seqkit`  - creates _.seqkit_ folder in your home directory.

`cd ~/.seqkit`   - change to .seqkit folder

`wget https://github.com/ashwini06/seqkit/blob/master/data/seqkit.yaml` - downloads the _seqkit.yaml_ file. 

Edit the root dir path in *seqkit.yaml*.
It should point to the path one level up where all the data folders are present

Eg:  `root_dir: "/proj/b2012025/RAW_DATA/ChIP_Data"`


To run seqkit we need some folder re-arrangements for placing the rawdata files. Seqkit considers [this](https://github.com/ashwini06/seqkit/blob/master/misc/Seqkit_project_dir.pdf) file structure. 

[Bash script](https://github.com/ashwini06/seqkit/blob/master/misc/arrange_data.sh) is written based on the general ChIP-Seq structure received during the time of project delivery. The script [rearranges](https://github.com/ashwini06/seqkit/blob/master/misc/SeqKit_projectdir_arrange.pdf) to the seqkit folder structure. However this script is not universal to all the project structures, as it is specifically written to one structure. 

**Note**
In the following example commands,  _Ascl1_US_ experiment is used to illustrate  [this](https://github.com/ashwini06/seqkit/blob/master/misc/Ascl1_US_projectstructure.pdf) structure. 

=======

### Available functions in seqkit
 
`seqkit  --help`

[Running fastqc](#preqc)

[Aligning and running bam-to-bed (default runs bowtie)](#align)

[Peak-calling (either in TF or HM mode: macs2/danpos2)](#Peak-call)

[Post-QC for analysis](#postqc)

[Runing STAR](#star)

[Running htSeq-Count, Cufflinks](#htcuff)

=======

<a name="preqc"></a>
### Running Fastqc 

A prior quality check to run fastqc on all the samples can be analyzed through _preqc_. It is important to check the sequence quality, over-represented sequences and duplicate percentages from the fastqc output results.
Fastqc output explaination is summarized [here](http://www4.ncsu.edu/~rosswhet/BIT815/Overview/Week2/FastQC_details.pdf)

`seqkit preqc -p Ascl1_US`


**Inputs**

| Parameters | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |
|-o/--out_dir | FILENAME (optional) | Path to output folder |


**Outputs**

Default : Creates fastqc folder inside the sample folder. _*.html_ contains the fastqc summarized report.
If _out_dir_ is specified _fastqc_ results will be stored in the desired path. 

=======

<a name="align"/></a>
### Aligning and running bam-to-bed (default runs bowtie)

[*Picard tools*](http://broadinstitute.github.io/picard/) is incorporated into _seqkit_ to estimate the quality metrics and [remove the duplicates](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) from aligned bam files.
However both the bam files _(*_sorted.bam)_ and _(*_rmdup.bam)_ files are present in the _(alignment_*)_ folder.
For further steps in seqkit, duplicates removed from the bam files are used .

**To run on all samples present in the project folder (-p)** 

`seqkit analysis align -p Ascl1_US --bam_to_bed`

Bed files are input to the peak calling software (macs2/danpos2).
So bed files can be generated while aligning the reads by adding an extra option *--bam_to_bed*.

Also, the above can be independently run (in case if we dont want to generate *bed files* 
or if we want to generate *bed files* from already existed *aligned reads*)

**Command-line**

`seqkit analysis align -p Ascl1_US ` (Doesnt generate bed files)

`seqkit analysis bam_to_bed -p Ascl1_US --slurm` (Generates bed files from already existed files)

**To run on specific samples (s) present in the project folder (-p)**

`seqkit analysis align -p Ascl1_US --bam_to_bed -s Mark_Mash1_s4`

`seqkit analysis align -p Ascl1_US --bam_to_bed -s Mark_input_s3`


**Inputs**

|  Parameters | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |
| -s/--sample | FILENAME *(optional)*| to run on specific samples inside project folder |
|-a/--aligner|bowtie/bowtie2/bwa |  which sofware to be used for aligning raw reads (default (bowtie))|
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

=======

<a name="Peak-call"/></a>

### Peak-calling (either in TF or HM mode: macs2/danpos2)

For _TF chip_, to call the enriched regions _macs2_ works well.
Whereas for some of the broad marks, macs2 might not work well 
(the called regions is divided into small segmental peaks).
Some of the publications [ref] used danpos2 to call the enriched regions for HM data.

**Command-line**

 `seqkit peakanalysis peakcall -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt`

## To annotate peak-called regions while peakcalling

`seqkit peakanalysis peakcall -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt --peakannotate`


**Inputs**

|  Parameters | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |
| -i/--input | FILENAME | Path to tab-delimited text file. It is a 2 column file, where first column should contain sample name (Treatment) and second column should contain input name (Control) |
|-m/--mode | TF/HM |  TF chip or Histone modification (default is TF) |
|--peak_call| macs2/danpos2|macs2 (applicable for TF and HM) or danpos2(only for HM) (default is HM) |

**Outputs**

[MACS outputs and parameters](https://github.com/taoliu/MACS)

[danpos2 outputs and parameters](https://sites.google.com/site/danposdoc/tutorial/dpeak)

**Note**

macs2/danpos2 command line parameters can be edited in [configuration file](https://github.com/ashwini06/seqkit/blob/master/data/seqkit.yaml)


## Seperate peak-annotations on already called peak-regions

` seqkit peakanalysis peakanno -p Ascl1_US --slurm`


=======

<a name="postqc"/></a>
### Post-QC for analysis

[*DeepTools*] (http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html) provides number of quality metrics and also an estimate for assessing the quality of ChIP. Incoporated _bamcompare_, _computeMatrix_, _plotHeatmap_, _plotCorrelation_, _plotFingerprint_ functions, which gives the idea about the chip enrichment signals, distribution of reads across genes, correlation between samples and whether ChIP experiment worked or not.

[*ngsplot*](https://github.com/shenlab-sinai/ngsplot) is another tool to create enrichment plots. This function complements with _plotHeatmap_ function from _DeepTools_ . 

**Command-line**

`seqkit postqc -p Ascl1_US -i /home/ashwini/scripts/seqkit/data/test.txt --genefile /home/ashwini/scripts/seqkit/data/UCSC_mm10genes.bed`

**Inputs**

|  Parameters | Expected Input | Explanation |
|:----:|:----:|:----|
| -p/--project | FILENAME	| Path to the project folder |
| -i/--input | FILENAME | Path to tab-delimited text file. It is a 2 column file, where first column should contain sample name (Treatment) and second column should contain input name (Control) |
|--genefile| FILENAME | Path to File name in BED format, containing the regions to plot |

**Outputs**

| Output files| Description|
|:----:|:----:|:----|
| *_coverage.bw | BigWig files can be uploaded to UCSC genome browser to view the coverage.Easiest way is to upload to the local server and use http link in the [track line](https://genome.ucsc.edu/goldenpath/help/bigWig.html). [Example](https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A15367354-15375823&hgsid=498762981_ascNmruLMMNBkmbedPsfE5wxlZiw)|
| *_fingerprint.png | Gives the overview of ChIPseq files. [More details](http://deeptools.readthedocs.io/en/latest/content/tools/plotFingerprint.html). [Example](https://github.com/ashwini06/seqkit/blob/master/data/Ascl1_Example/Mark_Mash1_s4_Vs_Mark_input_s3_fingerprint.png)|
| scatterplot.pdf   | Correlation value. [More details](http://deeptools.readthedocs.io/en/latest/content/tools/plotCorrelation.html).[Example](https://github.com/ashwini06/seqkit/blob/master/data/Ascl1_Example/scatterplot.pdf)|
|*avgprof.pdf| Plots the average profile across the genomic regions. Example [Average_Profile](https://github.com/ashwini06/seqkit/blob/master/data/Ascl1_Example/Mark_Mash1_s4_Vs_Mark_input_s3.genebody.avgprof.pdf)|
|*heatmap.pdf| Plots the heatmap for scores associated with genomic regions. [More details](http://deeptools.readthedocs.io/en/latest/content/tools/plotHeatmap.html).  Example  [heatmap](https://github.com/ashwini06/seqkit/blob/master/data/Ascl1_Example/Mark_Mash1_s4_Vs_Mark_input_s3.genebody.heatmap.pdf)|




## RNA-Seq analysis

<a name="star"></a>
### Running STAR

STAR is an ultra-fast aligner, which is generally used for aligning RNA-Seq reads

`seqkit analysis align -p lizzy_andersson_ChIP_ESC_57-64 -a STAR -s 130917_D3.5_wt_Rep1_H3K4me1`

<a name="htcuff"></a>
### Generating counts and fpkm values

[htSeq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) and [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html) are used to generate the counts and fpkm values. 

`seqkit htcuff -p lizzy_andersson_ChIP_ESC_57-64 -a STAR -s 130917_D3.5_wt_Rep1_H3K4me1`
