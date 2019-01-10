#!/bin/bash                                                                                                                                       

#SBATCH --job-name=MACS2_PeakCalling                                                             
##SBATCH --nodes=1                                                                                                              
##SBATCH --cpus-per-task=1                                                                                                                        
#SBATCH --mem=100GB                                                                                                                               
##SBATCH --gres=gpu:1                                                                                                                             
##SBATCH --partition=gpu4_medium                                                                                                                  
#SBATCH --partition=fn_short                                                                                                                      
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny/chipseq/20181217_FCHMYTCBGX7/err_out/%x_%j.err                                           
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/chipseq/20181217_FCHMYTCBGX7/err_out/%x_%j.out                                          
##SBATCH --dependency=afterany:job_id

module purge
module load macs2/2.1.1

# bam files
treatment=$1
ctr=$2
outdir=$3

IFS='/' read -ra ADDR <<< "$treatment"
treat_filename=${ADDR[-1]}
echo "treatment_filename is $treat_filename"
IFS='/' read -ra ADDR <<< "$ctr"
ctr_filename=${ADDR[-1]}
echo "control_filename is $ctr_filename"

macs2 callpeak -t ${treatment} -c ${ctr} -f BAM -g hs -n PC_${treat_filename}_vs_${ctr_filename} -B -q 0.01 --outdir ${outdir}
