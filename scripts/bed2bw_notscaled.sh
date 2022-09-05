#!/bin/bash
bed=${1?please give .bed file}
out=${2?output prefix}
genomesize=${3?please give chromosize file}

bedtools genomecov -i $bed -g $genomesize -bga | sort -k 1,1 -k 2,2n > ${out}.bdg
bedGraphToBigWig ${out}.bdg $genomesize ${out}.bw
rm ${out}.bdg

