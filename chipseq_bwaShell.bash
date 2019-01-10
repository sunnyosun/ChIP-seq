#!/bin/bash
#$ -S /bin/bash                                                                                                                                                                 
#$ -cwd

wd=`pwd`

for i in $wd/*.fastq; do qsub $wd/chipseq_bwa.bash $i; done