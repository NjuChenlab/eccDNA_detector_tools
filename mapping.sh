#!/bin/bash

outDir=./
thread=1
a1=AGATCGGAAGAGCACACGTC
a2=AGATCGGAAGAGCGTCGTGT
out=out

BASEDIR=$(dirname "$0")
script=$BASEDIR/scripts/

usage(){
  cat <<EOF
Usage: $(basename $0) Options
Options:
  -1|--r1           PATH    read 1 in fastq(.gz) format
  -2|--r2           PATH    read 2 in fastq(.gz) format
  -a|--adaptor      CHAR    3' adapter to be removed from read 1 in a pair [$a1]
  -A|--ADAPTOR      CHAR    3' adapter to be removed from read 2 in a pair [$a2]
  -i|--index        PREFIX  BWA Index prefix
  -o|--outDir       PATH    Output directory [$outDir]
  -p|--prefix       CHAR    prefix of output files [$out]
  -t|--thread       INT     The thread number [$thread]
  -h|--help
Dependencies:
  SeqPrep; samblaster v0.1.25; bwa v0.7.17-r1188; samtools v1.3.1
EOF
}

if [ $# -eq 0 ]; then
  usage; exit 1
fi

while [ -n "$1" ];do
  case $1 in
    -1|--r1)            shift;fq1=$1;;
    -2|--r2)            shift;fq2=$1;;
    -a|--adaptor)       shift;a1=$1;;
    -A|--ADAPTOR)       shift;a2=$1;;
    -i|--index)         shift;genome=$1;;
    -o|--outDir)        shift;outDir=$(readlink -mn $1);;
    -p|--prefix)        shift;out=$1;;
    -t|--thread)        shift;threads=$1;;
    -h|--help)          usage;exit 1;; 
    *)                  usage;exit 1;;
  esac
  shift
done

if [ -z $fq1 ];then
  echo "Error: fastq file for read 1 is not provided"; exit 1
fi
if [ ! -e $fq1 ];then
  echo "Error: fastq file for read 1: $fq1 doesn't exist"; exit 1
fi
if [ -z $fq2 ];then
  echo "Error: fastq file for read 2 is not provided"; exit 1
fi
if [ ! -e $fq2 ];then
  echo "Error: fastq file for read 2: $fq2 doesn't exist"; exit 1
fi
if [ -z $genome ];then
  echo "Error: prefix for BWA Index is not provided"; exit 1
fi
if [ ! -e ${genome}.bwt ];then
  echo "Error: prefix for BWA Index: $genome likely doesn't exist"; exit 1
fi

if ! command -v SeqPrep &> /dev/null
then
  echo "Error: SeqPrep could not be found"; exit 1
fi
if ! command -v samblaster &> /dev/null
then
  echo "Error: samblaster could not be found"; exit 1
fi
if ! command -v bwa &> /dev/null
then
  echo "Error: bwa could not be found"; exit 1
fi
if ! command -v samtools &> /dev/null
then
  echo "Error: samtools could not be found"; exit 1
fi


echo ">> preparing fastq files..."
echo ">> merging paired end Illumina reads that are overlapping into a single longer read..."
mkdir -p $outDir/map/${out}_fq
minlen=50
echo ">> only reads with length of >= $minlen are kept."
echo ">> log information for SeqPreq:"
SeqPrep -A $a1 -B $a2 -f $fq1 -r $fq2 \
        -1 $outDir/map/${out}_fq/${out}.r1.fq.gz -2 $outDir/map/${out}_fq/${out}.r2.fq.gz -L $minlen -s $outDir/map/${out}_fq/${out}.merge.fq.gz

echo ">> reads mapping by BWA..."
echo ">> reads mapping for paired end reads..."
bwa mem -q -M -v 1 -t $threads $genome $outDir/map/${out}_fq/${out}.r1.fq.gz $outDir/map/${out}_fq/${out}.r2.fq.gz | \
  samblaster -r -M -d $outDir/map/${out}.pe.disc.sam -s $outDir/map/${out}.pe.split.sam | \
  samtools view -Sb -o $outDir/map/${out}.pe.rmdup.bam
samtools view -bS $outDir/map/${out}.pe.disc.sam  > $outDir/map/${out}.pe.disc.bam  && rm $outDir/map/${out}.pe.disc.sam
samtools view -bS $outDir/map/${out}.pe.split.sam > $outDir/map/${out}.pe.split.bam && rm $outDir/map/${out}.pe.split.sam
echo ">> reads mapping for long singleton reads..."
bwa mem -q -M -v 1 -t $threads $genome $outDir/map/${out}_fq/${out}.merge.fq.gz | samtools view -Sb -o $outDir/map/${out}.se.bam
${script}/rmdup_sese.pl -i $outDir/map/${out}.se.bam | samtools view -Sb -o $outDir/map/${out}.se.rmdup.bam
echo ">> FINISHED! CHEERS!"
