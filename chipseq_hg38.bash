#!/bin/bash
#$ -S /bin/bash                                                                                                                                                                 
#$ -cwd

module load bwa
module load samtools
module load bedtools

fq_r1=$1.fastq.gz
sam=$1.sam
bam=$1.bam
#sortbam=$1.sorted.bam

echo "fq_r1=${fq_r1}"
echo "sortbam=$sortbam"

bwa mem -t 4 /ifs/data/proteomics/projects/Sunny/genome/hg38.fa /ifs/data/proteomics/projects/Sunny/fastq_chip_new/${fq_r1} | samtools view -b -S - > ${bam}

#samtools sort -O bam -T /tmp/$1 -m 2G -l 9 - -o ${sortbam}
#samtools index ${sortbam}

#bedtools genomecov -ibam ${sortbam} -bg > $1-coverage.bedgraph