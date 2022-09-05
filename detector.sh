#!/bin/bash

BASEDIR=$(dirname "$0")
overhang=10
quality=20
dist=500
script=$BASEDIR/scripts/
outDir=./
size=5000000

usage(){
    cat <<EOF
Usage: $(basename $0) Options
Options:
  -i|--input        BASE    prefix of the mapping.sh output
  -q|--mapq         INT     min mapping quality [20]
  -v|--overhang     INT     min length of overhangs beyond homology regions for lconf split reads detection [10]
  -g|--genome       FILE    genome sequence in fasta format
  -s|--size         FILE    text file for chromosome sizes
  -d|--distance     INT     max distance between 2 properly-paired / discordant reads [500]
  -l|--largest      INT     candidate eccDNA larger than [-l|--largest] will be disgarded [5000000]
  -o|--outDir       PATH    Output directory [$outDir]
  -h|--help
Dependencies:
  samtools; bedtools  
EOF
}

if [ $# -eq 0 ]; then
  usage; exit 1
fi

while [ -n "$1" ];do
  case $1 in
    -i|--input)         shift;input=$1;;
    -q|--mapq)          shift;quality=$1;;
    -v|--overhang)      shift;overhang=$1;;
    -g|--genome)        shift;genome=$1;;
    -s|--size)          shift;chromsizes=$1;;
    -d|--distance)      shift;dist=$1;;
    -l|--largest)       shift;size=$1;;
    -o|--outDir)        shift;outDir=$(readlink -mn $1);;
    -h|--help)          usage;exit 1;; 
    *)                  usage;exit 1;;
  esac
  shift
done

samp=`basename $input`

if [ -z $input ];then
  echo "Error: -i|--input is not provided"; exit 1
fi
if [ ! -e ${input}.pe.split.bam ];then
  echo "Error: ${samp}*.bam files likely do not exist"; exit 1
fi
if [ -z $genome ];then
  echo "Error: genome sequence in .fa format is not provided"; exit 1
fi
if [ ! -e $genome ];then
  echo "Error: $genome doesn't exist"; exit 1
fi
if [ -z $chromsizes ];then
  echo "Error: file for chromosome sizes is not provided"; exit 1
fi
if [ ! -e $chromsizes ];then
  echo "Error: $chromsizes doesn't exist"; exit 1
fi
if ! command -v bedtools &> /dev/null
then
  echo "Error: bedtools could not be found"; exit 1
fi
if ! command -v samtools &> /dev/null
then
  echo "Error: samtools could not be found"; exit 1
fi

echo ">> detecting hconf split reads for paired-end reads..."
echo "  >> detecting candidate hconf split reads for r1 or r2 separatedly..." 
mkdir -p $outDir/out
$script/detect_hconf_split.pl -i ${input}.pe.split.bam -q $quality | \
  awk -v l=$size '$8-$5<=l' > $outDir/out/${samp}.pe.out

echo "    >> both r1 and r2 are hconf split reads..." 
# homology sequences around break points are taken into consideration
# ........ATCGATCG|................................ATCGATCG|.........
# OR
# .......|ATCGATCG................................|ATCGATCG..........
#              1              2 3  4      5        6        7        8     9              0              1 2   3      4        5        6       7     8
# chrX:10901490-10901609:15-1 1 + chrX 10901491 10901556 10901575 10901610 1 chrX:10901490-10901609:15-1 2 - chrX 10901491 10901539 10901558 10901610 1
cut -f 1 $outDir/out/${samp}.pe.out | sort | uniq -d | \
  $script/filter.pl -o - -1 1 -2 1 -m i $outDir/out/${samp}.pe.out | sed 'N;s/\n/\t/' | \
  awk '$3!=$12 && $4==$13' | awk '($5==$14 && $8==$17)' | \
  awk 'BEGIN{OFS="\t"}{print $1,12,$4,$5,max($6,$15),min($7,$16),$8,$9}
       function max(a, b){return a > b ? a : b}
       function min(a, b){return a > b ? b : a}' > $outDir/out/${samp}.pe.flt.out
cut -f 1 $outDir/out/${samp}.pe.out | sort | uniq -d | \
  $script/filter.pl -o - -1 1 -2 1 -m i $outDir/out/${samp}.pe.out | sed 'N;s/\n/\t/' | \
  awk '$3!=$12 && $4==$13' | awk '($5!=$14 || $8!=$17)' | awk '$9>$18' | \
  awk '$14+$9-$18==$5 && $17+$9-$18==$8' | \
  awk 'BEGIN{OFS="\t"}{print $1,12,$4,$5,max($6,$15),min($7,$16),$8,$9}
       function max(a, b){return a > b ? a : b}
       function min(a, b){return a > b ? b : a}' >> $outDir/out/${samp}.pe.flt.out
cut -f 1 $outDir/out/${samp}.pe.out | sort | uniq -d | \
  $script/filter.pl -o - -1 1 -2 1 -m i $outDir/out/${samp}.pe.out | sed 'N;s/\n/\t/' | \
  awk '$3!=$12 && $4==$13' | awk '($5!=$14 || $8!=$17)' | awk '$9<$18' | \
  awk '$5+$18-$9==$14 && $8+$18-$9==$17' | \
  awk 'BEGIN{OFS="\t"}{print $1,12,$4,$14,max($6,$15),min($7,$16),$17,$9}
       function max(a, b){return a > b ? a : b}
       function min(a, b){return a > b ? b : a}' >> $outDir/out/${samp}.pe.flt.out

echo "    >> r1 (or r2) is a split read, and r2 (or r1) is non-split read"
# get the other read
cut -f 1 $outDir/out/${samp}.pe.out | sort | uniq -u | \
  $script/filter.pl -o - -1 1 -2 1 -m i <(samtools view ${input}.pe.rmdup.bam) | \
  cat <(samtools view -H ${input}.pe.rmdup.bam) - | samtools view -Sb | \
  bedtools bamtobed | sed 's/\//\t/g' | awk -v qua=$quality '$6>=qua' | \
  $script/filter.pl -o $outDir/out/${samp}.pe.out -1 1,2 -2 4,5 -s : -m e > $outDir/out/${samp}.pairedread
# uniquely-mapped;
# mapped on the opposite strand;
# mapped on the same chromosome;
# within the eccDNA body (specially, 20nt-extra-region is added to the eccDNA to accommodate additional homology sequences flanking the eccDNA)
# dist <= $dist
cut -f 4 $outDir/out/${samp}.pairedread | sort | uniq -d | $script/filter.pl -o - -1 1 -2 4 -m e $outDir/out/${samp}.pairedread | \
  $script/skyjoin -1 1 -2 4 $outDir/out/${samp}.pe.out - | awk '$3!=$15' | awk '$4==$10' | awk '($5-$9 - 20) <=$11 && ($8 + 20)>=$12' | \
  awk -v d=$dist '($3=="+" && $12-$5+$8-$7<=d) || ($3=="-" && $6-$5+$8-$11<=d)' | \
  awk 'BEGIN{OFS="\t"}{print $1,12,$4,$5,$6,$7,$8,$9}' >> $outDir/out/${samp}.pe.flt.out
rm $outDir/out/${samp}.pairedread

echo ">> detecting hconf split reads for longer singleton reads..."
samtools view ${input}.se.rmdup.bam | cut -f 1 | uniq -c | sed 's/^  *//g' | sed 's/  */\t/g' | awk '$1==2' | \
  $script/filter.pl -o - -1 2 -2 1 -m i <(samtools view ${input}.se.rmdup.bam) | \
  $script/detect_hconf_split.pl -q $quality | awk -v l=$size '$8-$5<=l' | \
  awk 'BEGIN{OFS="\t"}{print $1,12,$4,$5,$6,$7,$8,$9}' > $outDir/out/${samp}.se.flt.out

echo ">> combining split reads detected from paired-ends reads and singleton reads..."
# all individual read pairs with break points detected (left_least_start, right_most_end)
cat $outDir/out/${samp}.pe.flt.out $outDir/out/${samp}.se.flt.out > $outDir/out/${samp}.bks.indi
# reads ids examined
cat $outDir/out/${samp}.bks.indi | cut -f 1 | sort -u > $outDir/out/${samp}.id.passed
rm  $outDir/out/${samp}.pe.flt.out
rm  $outDir/out/${samp}.se.flt.out
# lconf
$script/filter.pl -o $outDir/out/${samp}.bks.indi -1 1 -2 1 -m e $outDir/out/${samp}.pe.out > $outDir/out/${samp}.lconf.out
rm  $outDir/out/${samp}.pe.out

echo ">> reporting candidate eccDNAs supported by hconf split reads..."
echo "  >> format: chr left_most_start right_least_end name offset num_of_hconf_split_reads names_of_hconf_split_reads"
$script/skycut.pl -f 3,4,7,1,8 $outDir/out/${samp}.bks.indi | sort -k 1,1 -k 2,2n -k 3,3n -k 5,5n | bedtools groupby -g 1,2,3,5 -c 4,4 -o count,collapse | \
  awk 'BEGIN{OFS="\t"}{print $1,$2-$4,$3-$4,$4,$5,$6}' | \
  sort -k 1,1 -k 2,2n -k 3,3n | bedtools groupby -g 1,2,3 -c 4,5,6 -o max,sum,collapse | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,NR,$4,$5,$6}'> $outDir/out/${samp}.candidate.eccDNA

echo ">> retriving lconf split reads..."
echo "  >> for properly-paired read pairs..."
# don't need to re-check those split reads
$script/filter.pl -o $outDir/out/${samp}.id.passed -1 1 -2 1 -m e <(samtools view -f 0x2 -F 0x100 -F 0x800 ${input}.pe.rmdup.bam) | \
  $script/detect_lconf_split_con.pl -e $outDir/out/${samp}.candidate.eccDNA -q $quality -o $overhang -f $genome | \
  sort -k 1,1n > $outDir/out/${samp}.pe.con.lconf
echo "  >> for discordant read pairs..."
$script/filter.pl -o $outDir/out/${samp}.id.passed -1 1 -2 1 -m e <(samtools view -F 0x100 -F 0x800 ${input}.pe.disc.bam) | \
  $script/detect_lconf_split_disc.pl -e $outDir/out/${samp}.candidate.eccDNA -q $quality -o $overhang -f $genome -d $dist | \
  sort -k 1,1n > $outDir/out/${samp}.pe.disc.lconf
echo "  >> for singleton reads..."
$script/filter.pl -o $outDir/out/${samp}.id.passed -1 1 -2 1 -m e <(samtools view -F 0x100 -F 0x800 ${input}.se.rmdup.bam) | \
  $script/detect_lconf_split_se.pl -e $outDir/out/${samp}.candidate.eccDNA -q $quality -o $overhang -f $genome | \
  sort -k 1,1n > $outDir/out/${samp}.se.lconf

echo ">> reporting lconf split reads for individual candidate eccDNAs..."
paste $outDir/out/${samp}.pe.con.lconf $outDir/out/${samp}.pe.disc.lconf $outDir/out/${samp}.se.lconf | \
  awk 'BEGIN{OFS="\t"}{print $1,$2+$5+$8,$3,$6,$9}' | \
  sed 's/\t/,/3' | sed 's/\t/,/3' | sed 's/-,//g' | sed 's/-$//g' | sed 's/,$//g' | sed 's/\t$/\t-/g' | \
  paste $outDir/out/${samp}.candidate.eccDNA - | \
  $script/skycut.pl -f 1,2,3,4,5,6,9,7,10 > $outDir/out/${samp}.candidate.eccDNA.full
mv $outDir/out/${samp}.candidate.eccDNA.full $outDir/out/${samp}.candidate.eccDNA
cut -f 9 $outDir/out/${samp}.candidate.eccDNA | grep -v -P "^-$" | sed 's/,/\n/g' >> $outDir/out/${samp}.id.passed
rm $outDir/out/${samp}.pe.con.lconf
rm $outDir/out/${samp}.pe.disc.lconf
rm $outDir/out/${samp}.se.lconf

echo ">> detecting supporting disc reads for individual candidate eccDNAs..."
$script/filter.pl -o $outDir/out/${samp}.id.passed -1 1 -2 1 -m e <(samtools view -F 0x100 -F 0x800 ${input}.pe.disc.bam) | \
  cat <(samtools view -H ${input}.pe.disc.bam) - | samtools view -Sb | bedtools bamtobed -cigar | sed 'N;s/\n/\t/' | awk '$1==$8 && $6!=$13' | \
  awk 'BEGIN{OFS="\t"}($6=="+" && ($2>$9 || $3>$10)) || ($6=="-" && ($2<$9 || $3<$10)){print $0}' | \
  awk -v qual=$quality '$5>=qual && $12>=qual' | \
  bedtools intersect -a <(cut -f 1-7 $outDir/out/${samp}.candidate.eccDNA) -b - -wo | \
  awk -v dist=$dist 'BEGIN{OFS="\t"}min($9,$16)>=$2 && max($10,$17)<=$3 && (min($9,$16)-$2+$3-max($10,$17))<=dist{print $1,$2,$3,$4,$5,$6,$7,$11}
       function min(a, b){return a > b ? b : a}function max(a, b){return a > b ? a : b}' | \
  sed 's/\/1//g' |  sort | bedtools groupby -g 1,2,3,4,5,6,7 -c 8,8 -o count_distinct,distinct > $outDir/out/${samp}.candidate.eccDNA.wdisc

$script/filter.pl -o $outDir/out/${samp}.candidate.eccDNA.wdisc -1 4 -2 4 -m e $outDir//out/${samp}.candidate.eccDNA | \
  cut -f 1-7 | sed 's/$/\t0\t-/g' | cat - <(cut -f 1-9 $outDir/out/${samp}.candidate.eccDNA.wdisc) | sort -k 4,4n | \
  cut -f 8,9 | paste <(sort -k 4,4n $outDir/out/${samp}.candidate.eccDNA) - | \
  $script/skycut.pl -f 1-7,10,8,9,11 > $outDir/out/${samp}.eccDNA+reads
rm $outDir/out/${samp}.candidate.eccDNA $outDir/out/${samp}.candidate.eccDNA.wdisc

echo ">> post-detection steps..."

echo "  >> preparing .bed files for split reads..."
echo "    >> for hconf split reads"
cat $outDir/out/${samp}.eccDNA+reads | cut -f 9 | sed 's/,/\n/g' | sort -u | \
  $script/filter.pl -o - -1 1 -2 1 -m i <(cat <(samtools view ${input}.pe.rmdup.bam) <(samtools view ${input}.se.rmdup.bam)) | \
  cat <(samtools view -H ${input}.pe.rmdup.bam) - | samtools view -Sb | \
  bedtools bamtobed -cigar | sed 's/\/[12]//g' | awk -v qual=$quality '$5>=$qual' > $outDir/out/${samp}.split.bed
echo "    >> for lconf split reads"
cat $outDir/out/${samp}.eccDNA+reads | cut -f 10 | sed 's/,/\n/g' | sed '/^-$/d' | sort -u | \
  $script/filter.pl -o - -1 1 -2 1 -m i <(cat <(samtools view -F 0x100 -F 0x800 ${input}.pe.rmdup.bam) <(samtools view -F 0x100 -F 0x800 ${input}.se.rmdup.bam)) | \
  cat <(samtools view -H ${input}.pe.rmdup.bam) - | samtools view -Sb | \
  bedtools bamtobed -cigar | sed 's/\/[12]//g' | awk -v qual=$quality '$5>=$qual' >> $outDir/out/${samp}.split.bed

echo "  >> preparing .bed files for disc reads..."
cut -f 11 $outDir/out/${samp}.eccDNA+reads | sed 's/,/\n/g' | sort -u | \
  $script/filter.pl -o - -1 1 -2 1 -m i <(samtools view ${input}.pe.disc.bam) | \
  cat <(samtools view -H ${input}.pe.disc.bam) - | samtools view -Sb | \
  bedtools bamtobed -cigar | sed 's/\/[12]//g' > $outDir/out/${samp}.disc.bed

echo "  >> preparing .bed files for properly-paired non-split reads..."
cut -f 9,10,11 $outDir/out/${samp}.eccDNA+reads | tr "\t" "\n" | sed 's/,/\n/g' | sed '/^-$/d' | \
  $script/filter.pl -o - -1 1 -2 1 -m e <(samtools view -h -f 0x2 -F 0x100 -F 0x800 ${input}.pe.rmdup.bam) | \
  samtools view -Sb | bedtools bamtobed -cigar | sed 'N;s/\n/\t/' | awk -v qual=$quality '$5>=qual && $12>=qual' | \
  sed 's/\t/\n/7' | sed 's/\/[12]//g' > $outDir/out/${samp}.con.bed
cut -f 9,10,11 $outDir/out/${samp}.eccDNA+reads | tr "\t" "\n" | sed 's/,/\n/g' | sed '/^-$/d' | \
  $script/filter.pl -o - -1 1 -2 1 -m e <(samtools view -h -F 0x100 -F 0x800 ${input}.se.rmdup.bam) | \
  samtools view -Sb | bedtools bamtobed -cigar | awk -v qual=$quality '$5>=qual' | \
  sed 's/\/[12]//g' >> $outDir/out/${samp}.con.bed
echo "  >> combining all .bed files"
sort -k 1,1 -k 2,2n -k 3,3n $outDir/out/${samp}.split.bed $outDir/out/${samp}.disc.bed $outDir/out/${samp}.con.bed > $outDir/out/${samp}.bed
rm $outDir/out/${samp}.split.bed $outDir/out/${samp}.disc.bed $outDir/out/${samp}.con.bed

echo "  >> reporting all covered regions..."
bedtools merge -i $outDir/out/${samp}.bed > $outDir/out/${samp}.cov

echo "  >> calculating %covered for individual eccDNA and average coverage in RPK..."
cut -f 1-8 $outDir/out/${samp}.eccDNA+reads | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3+$5,$4,$5,$6,$7,$8}' | \
  bedtools intersect -a - -b $outDir/out/${samp}.cov -wo | cut -f 1-8,12 | \
  bedtools groupby -g 1-8 -c 9 -o sum | awk 'BEGIN{OFS="\t"}{print $0,$9/($3-$2)}' | cut -f 1-8,10 | \
  bedtools intersect -a $outDir/out/${samp}.bed -b - -wo | $script/skycut.pl -f 8-16,4 | sort | \
  bedtools groupby -g 1-9 -c 10 -o count_distinct | awk 'BEGIN{OFS="\t"}{print $0,$10/($3-$2)*1000}' | \
  cut -f 1-9,11 | sort -k 4,4n | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3-$5,$4,$5,$6,$7,$8,$9,$10}' > $outDir/out/${samp}.eccDNA

echo "  >> calculating average coverage of upstream regions..."
# add a few bp offset (e.g., 10bp) to more separate with eccDNA
cut -f 1-4 $outDir/out/${samp}.eccDNA+reads | \
  bedtools slop -i - -g $chromsizes -l 10 -r -10 | \
  bedtools flank -i - -g $chromsizes -l 1 -r 0 -pct | cut -f 1-4 | awk '$3-$2!=0' > $outDir/out/${samp}.l.bed
bedtools intersect -a $outDir/out/${samp}.bed -b $outDir/out/${samp}.l.bed -wo | $script/skycut.pl -f 8-11,4 | sort | \
  bedtools groupby -g 1-4 -c 5 -o count_distinct | sort -k 4,4n > $outDir/out/${samp}.l
$script/filter.pl -o $outDir/out/${samp}.l -1 4 -2 4 -m e $outDir/out/${samp}.l.bed | sed 's/$/\t0/g' | cat $outDir/out/${samp}.l - | \
  awk 'BEGIN{OFS="\t"}{print $4,$5/($3-$2)*1000}' | \
  cat - <($script/filter.pl -o $outDir/out/${samp}.l.bed -1 4 -2 4 -m e $outDir/out/${samp}.eccDNA | cut -f 4 | sed 's/$/\tna/g') | \
  sort -k 1,1n > $outDir/out/${samp}.l.temp
mv $outDir/out/${samp}.l.temp $outDir/out/${samp}.l 
rm $outDir/out/${samp}.l.bed 

echo "  >> calculating average coverage of downstream regions..."
# add a few bp offset (e.g., 10bp) to more separate with eccDNA
cut -f 1-5 $outDir/out/${samp}.eccDNA+reads | \
  awk 'BEGIN{OFS="\t"}{print $1,$2+$5,$3+$5,$4}' | \
  bedtools slop -i - -g $chromsizes -l -10 -r 10 | \
  bedtools flank -i - -g $chromsizes -l 0 -r 1 -pct | cut -f 1-4 | awk '$3-$2!=0' > $outDir/out/${samp}.r.bed
bedtools intersect -a $outDir/out/${samp}.bed -b $outDir/out/${samp}.r.bed -wo | $script/skycut.pl -f 8-11,4 | sort | \
  bedtools groupby -g 1-4 -c 5 -o count_distinct | sort -k 4,4n > $outDir/out/${samp}.r
$script/filter.pl -o $outDir/out/${samp}.r -1 4 -2 4 -m e $outDir/out/${samp}.r.bed | sed 's/$/\t0/g' | cat $outDir/out/${samp}.r - | \
  awk 'BEGIN{OFS="\t"}{print $4,$5/($3-$2)*1000}' | \
  cat - <($script/filter.pl -o $outDir/out/${samp}.r.bed -1 4 -2 4 -m e $outDir/out/${samp}.eccDNA | cut -f 4 | sed 's/$/\tna/g') | \
  sort -k 1,1n > $outDir/out/${samp}.r.temp
mv $outDir/out/${samp}.r.temp $outDir/out/${samp}.r 
rm $outDir/out/${samp}.r.bed

echo ">> final output format:"
echo ">> chr start end id offset hconf_split lconf_split disc coverage(0-1) cov(RPK) cov_upstream(RPK) cov_downstream(RPK)"
paste $outDir/out/${samp}.eccDNA $outDir/out/${samp}.l $outDir/out/${samp}.r | \
  cut -f 1-10,12,14 > $outDir/out/${samp}.eccDNA.temp
mv $outDir/out/${samp}.eccDNA.temp $outDir/out/${samp}.eccDNA
rm $outDir/out/${samp}.l $outDir/out/${samp}.r
echo ">> FINISHED! CHEERS!"
