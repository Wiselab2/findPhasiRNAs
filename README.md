# findPhasiRNAs

Detect genomic phased loci from small RNA data 

## AUTHOR/SUPPORT

Sagnik Banerjee, !(sagnik@iastate.edu) !(https://github.com/sagnikbanerjee15/findPhasiRNAs/issues)

## HARDWARE/SOFTWARE REQUIREMENTS

x86-64 compatible processors
64 bit Linux or Mac OS X

## INTRODUCTION

Plant genomes encode abundant but diverse populations of small non-coding RNAs, which can be broadly divided into microRNAs (miRNAs) and endogenous small interfering RNAs (siRNAs). Endogenous siRNAs can be further grouped into several sub-classes such as heterochromatic small interfering RNAs (siRNAs), natural antisense siRNAs (nat-siRNAs) and trans-acting siRNAs (ta-siRNAs). The role of miRNAs as post-transcriptional regulators is well known. Among siRNAs, tasiRNAs and natsiRNAs are known to act as guide molecules for post-transcriptional gene regulation, and heterochromatic siRNAs in transcriptional gene silencing, but the role of phasiRNAs in gene regulation is still unclear. PhasiRNAs are produced from both protein-coding and noncoding genes. In many eudicots, three large gene families generate the majority of phasiRNAs, including those encoding nucleotide binding leucine-rich repeat proteins (NB-LRR genes), pentatricopeptide repeat proteins (PPR genes), and MYB transcription factors (MYB genes). 

This pipeline takes a statistical approach to obtain genomic loci where there is a strong signal of phasing. It computes p-values and also phasing scores. Phasing structure for each potential phasic loci is output as an image file to the output directory.   

## How to run findPhasiRNAs

findPhasiRNAs is written mainly in Python. There is one R code to generate plots. Binaries of all other dependencies are included in the release. Please note that you must have python3 to run the code. Version 2 of python is not supported.

### Installation

**Install Python [Please skip of you already have python installed]**

Several distrbutions of Linux have python3 preinstalled. Please check have installed, open a command prompt and run

```
python3 --version
```

For Ubuntu 16.10 or Newer

```
sudo apt-get update
sudo apt-get install python3.6
sudo pip install biopython
```

If you’re using another version of Ubuntu (e.g. the latest LTS release)

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get install python3.6
sudo pip install biopython
```

Fedora or RHEL distributions

```
sudo dnf install python3
sudo pip install biopython
```

**Install R [Please skip if you already have R installed]**

For Ubuntu

```
sudo apt-get install r-base
```

Fedora or RHEL distributions

```
sudo yum install r
```

Install other dependencies using the following commands

```
install.packages("plyr")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("grid")
install.packages("gridExtra")
```

In case you are using a **SPACK** based cluster, you will need to load each module individually. Please note that the default versions could differ from what is provided below.

```
module load python
module load py-biopython/1.70-py3-wos466g
module load r
module load r-plyr
module load r-ggplot2
module load r-reshape2
module load r-gridextra
```

### Installing findPhasiRNAs

Apart from python and R no other softwares are required.

```
git clone https://github.com/sagnikbanerjee15/findPhasiRNAs.git
cd findPhasiRNAs
```
 
### Parameter description

The entire software has been written in python and has been put in a single file `findPhasiRNAs.py`. There is another file `plot.R` which plots the phasing score graphs.

Please trim adapters from your sequences. You could either use Trimmomatic !(http://www.usadellab.org/cms/?page=trimmomatic]) or cutadapt !(https://cutadapt.readthedocs.io/en/stable/guide.html).

The release file  
The program will require 3 mandatory inputs:
- Either the genome sequence or bowtie index of the genome sequence. 
- Input fasta file or consolidated counts file. You can provide the fastq file if you do not have the fasta file. The program will perform the required conversions.
- Output directory name. The program will create the output directory if it does not exist

Bowtie Index generation from the genome takes a long time. It is recommended to build the index before phasiRNA analysis if you have multiple samples. Constructing the indices once will considerably save time by eliminating redundant executions.  

```
python findPhasiRNAs.py --help
usage: findPhasiRNAs.py [-h]
                        (--input_library INPUT_LIBRARY | --consolidated_library CONSOLIDATED_LIBRARY)
                        (--genome GENOME | --bowtie_index BOWTIE_INDEX)
                        --output_directory_provided OUTPUT_DIRECTORY_PROVIDED
                        --small_rna_size SMALL_RNA_SIZE [SMALL_RNA_SIZE ...]
                        --number_of_cycles NUMBER_OF_CYCLES
                        [NUMBER_OF_CYCLES ...] --pvalue_cutoff PVALUE_CUTOFF
                        [--clean_up CLEAN_UP] [--CPU CPU]
                        [--map_limit MAP_LIMIT] [--force FORCE]

findPhasiRNAs can be used to find locations where phasing occurs. We recommend
that you trim adapters from your libraries before submitting them to this
pipeline. The pipeline will NOT perform any adapter trimming.

optional arguments:
  -h, --help            show this help message and exit
  --input_library INPUT_LIBRARY, -i INPUT_LIBRARY
                        Specify the name of the file which has the small-RNA
                        reads. This option is mutually exclusive with
                        --consolidated_library
  --consolidated_library CONSOLIDATED_LIBRARY, -clib CONSOLIDATED_LIBRARY
                        Specify the name of the file which has the reads
                        consolidated by the number of occurances. This must be
                        in fasta format. The fasta header of each read must be
                        followed by the number of times they occur in the
                        original dataset separated by an underscore. For
                        example, if the read occurs 90182 times, then the
                        fasta header should read <read_name>_90182. You can
                        provide any <read_name> as you wish. Please note that
                        the number of occurances of the reads will be used to
                        calculate phasing score. This option is mutually
                        exclusive with --input_library.
  --genome GENOME, -g GENOME
                        Specify the name of the genome fasta file of the
                        organism. Please note that the program will not be
                        able to handle multiple fasta files.
  --bowtie_index BOWTIE_INDEX, -bindex BOWTIE_INDEX
                        Provide the bowtie index. This argument is optional.
                        If no index is provided then the software will
                        generate one.

Optional Arguments:
  --small_rna_size SMALL_RNA_SIZE [SMALL_RNA_SIZE ...], -srnasize SMALL_RNA_SIZE [SMALL_RNA_SIZE ...]
                        Specify the size of the small RNA that you wish to
                        analyze. You can enter more than one possible size.
  --number_of_cycles NUMBER_OF_CYCLES [NUMBER_OF_CYCLES ...], -numcycles NUMBER_OF_CYCLES [NUMBER_OF_CYCLES ...]
                        Specify the number of cycles you wish to analyze with.
                        You can enter multiple number of number of cycles. The
                        accepted values are 9, 10, 11, 12 and 13
  --pvalue_cutoff PVALUE_CUTOFF, -p PVALUE_CUTOFF
                        Enter the p-value cut off
  --clean_up CLEAN_UP, -c CLEAN_UP
                        Set this to 1 if you wish to clean up all the
                        intermediate files. The program will keep all
                        temporary files by default.
  --CPU CPU, -n CPU     Provide the number of CPUs to be used. Default is 1.
  --map_limit MAP_LIMIT, -mapl MAP_LIMIT
                        Specify the mapping limit. Only reads which are mapped
                        at most -mapl times will be considered. The default is
                        1. The maximum number of alignments allowed for a
                        single read is 10.
  --force FORCE, -f FORCE
                        Overwrite contents of output directory if it exists.

Required Arguments:
  --output_directory_provided OUTPUT_DIRECTORY_PROVIDED, -out OUTPUT_DIRECTORY_PROVIDED
                        Specify an output directory to which all the generated
                        files will be housed. This includes the log file which
                        can be later checked. Please make sure that there are
                        sufficient permissions to create the output directory.
                        The program will throw an error if creation of the
                        output directory fails. If the directory already
                        exists then its contents will be overwritten without
                        warning. This directory will contain the summary file
                        containing the details of the execution
```

### Example Runs

Follow the following steps to run analysis with an example dataset

```
# Download the Arabidopsis Thaliana genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-42/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Download the fastq data from NCBI
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR510/SRR5100580/SRR5100580.sra
fastq-dump SRR5100580.sra

# Run with default arguments [Run1]
python findPhasiRNAs.py -i SRR5100580.fastq -g Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -out Run1  

```

## BACKGROUND FORMULA

### Algorithm for computation of p-value

* Step 1: Map sRNAs to reference genome using bowtie1. Perform exact match and remember to get rid of all information about those reads which do not map to any coordinate. This needs to be done to keep the alignment file size in check

* Step 2: Add a two-nucleotide positive offset to all the sRNAs mapped on to the anti-sense strand to account for the 3’ overhang of tasiRNAs. Combine two reads which appear on both strands in the same location.

* Step 3: Now we define some terms:

L: length of the sRNA you are interested in. This value is typically between 20 and 22.

_Phase Register_: A region of the genome where an sRNA could cleave. In some papers, this same thing is called a cycle, a phase cycle or a register. Typically, a phase register is as long as ‘L’.

_Window_: A sequence of the genome which will be probed for presence of sRNAs. This has a typical length of 9L, 10L or 11L. Length of a window (mL) is often represented as a multiple of the ‘m’ – the number of phase registers it contains.

![Diagram to explain sRNA p-value computation](images/Picture1.png)

The above figure depicts a window which in this case is 231bp long. Addition of 2 nt to coordinates of the antisense strand will align the lower strand end-to-end with the strand above it. Here L=21 and m=11.

Positive Windows: Windows that abide by the 3 following rules are called positive windows:
•	Contains at least 10 unique reads
•	More than half of the reads should be ‘L’ nt long
•	At least three unique reads falling into the phase registers (Not sure how important this point is – may choose to ignore it)

Phased and non-phased locations: The vertical arrow indicates the start site for the small RNA used to determine the phased and non-phased positions. 21 phased sites relative to the start site are indicated as black vertical bars. Four hundred forty non-phased sites relative to the start site are indicated as grey. In this paper [4], they have considered the two strands separately which is why there are more phased sites. In our case there will be ‘m’ phased sites in a window.

‘n’: Total number of possible locations where phasing can occur. Hence, in our case, a window can have maximum ‘m’ such positions.

‘k’: Number of phased locations in the window which is covered by at least one sRNA.

$$pvalue = \sum_{x=k}^{m}\frac{\binom{20m}{n-x}\binom{m}{x}}{\binom{21m}{n}}$$

phased sRNA: reads which have mapped onto one of the phased locations.

Step 4: Calculate the Phasing score

Phasing score will be computed for each location (loc) of the genome using the following formula.

$$PhaseScore_{loc} = ln(1+10 \times \frac{\sum_{i=1}^{m} P_i}{1+\sum_{i=1}^{m}U_i})^{k-2}$$

‘k’: Number of phased locations in the window which is covered by at least one sRNA. Or this calculation we will consider k>=3.
Pi: Number of phased reads at the ith phase from the position loc
Ui: total number of reads for all small RNAs with start coordinates out of the ith phase

## LIMITATIONS


## FUNDING


