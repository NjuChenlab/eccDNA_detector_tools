# eccDNA_detector_tools  
![advantage of this method](image/advantage.png)
## Introduction  
## Installation  
In order to utilize the eccDNA detector tool, you must first download the package listed in the table. We recommend using conda to install.  
Following list was the detailed dependence:
software  |version
|:---:  |:---:|
|SeqPrep|v. 1.2  |
| bwa | v. 0.7.17-r1188 |
|samblaster  |v. 0.1.26  |
|samtools  | v. 1.7 |
| bedtools |v. 2.30.0  |  

Then, you need to add the tool to your path.  
```
export PATH=$PATH:/home/username/eccDNA_detector_tools/  
source ~/.bashrc  
```
## Usage  
### step1 mapping  
```
mapping.sh -1 read1.fq.gz -2 read2.fq.gz  -i bwa_genome_index -o out_dir -p output_prefix -t 12  
Usage: mapping.sh Options
Options:
  -1|--r1           PATH    read 1 in fastq(.gz) format
  -2|--r2           PATH    read 2 in fastq(.gz) format
  -a|--adaptor      CHAR    3' adapter to be removed from read 1 in a pair [AGATCGGAAGAGCACACGTC]
  -A|--ADAPTOR      CHAR    3' adapter to be removed from read 2 in a pair [AGATCGGAAGAGCGTCGTGT]
  -i|--index        PREFIX  BWA Index prefix
  -o|--outDir       PATH    Output directory [./]
  -p|--prefix       CHAR    prefix of output files [out]
  -t|--thread       INT     The thread number [1]
  -h|--help
Dependencies:
  SeqPrep; samblaster v0.1.25; bwa v0.7.17-r1188; samtools v1.3.1

```
The output dir contain four type files
```
output_prefix.pe.disc.bam  output_prefix.pe.rmdup.bam  output_prefix.pe.split.bam  output_prefix.se.bam  output_prefix.se.rmdup.bam
```

### step2  Detecting  
```
detector.sh -i output_prefix -g genome.fa -s genome.size -o outputdir  
Usage: detector.sh Options
Options:
  -i|--input        BASE    prefix of the mapping.sh output
  -q|--mapq         INT     min mapping quality [20]
  -v|--overhang     INT     min length of overhangs beyond homology regions for lconf split reads detection [10]
  -g|--genome       FILE    genome sequence in fasta format
  -s|--size         FILE    text file for chromosome sizes
  -d|--distance     INT     max distance between 2 properly-paired / discordant reads [500]
  -l|--largest      INT     candidate eccDNA larger than [-l|--largest] will be disgarded [5000000]
  -o|--outDir       PATH    Output directory [./]
  -h|--help
Dependencies:
  samtools; bedtools
``
### step3
