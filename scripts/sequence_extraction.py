"""
This script will extract sequences from human genome.fa file based on start and end positions in a .bed file
by Xiaoji Sun
updated on 2/12/2019

Inputs:
1) genome file (e.g. hg38.fa)
2) query file (a bed file that specifies chr, start, end)
3) Optional: strand (index number of the column that contains strandness(+/-) infomation, if not specified, the sequence will be 
				extracted based on start and end. if start > end, reverse complement sequence will be exported.)
4) Optional: length (the minimum length of output sequences)
"""

# Usage: sequence_extraction.py [options]

################################################################################
# Modules

# import regular expression module
import re

# import sys module
import sys
import subprocess
import optparse
import logging
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd

logging.basicConfig(level=logging.INFO,
					format='%(asctime)s %(levelname)s %(message)s')

################################################################################
# Functions
# load in hg19 or hg38 and extract only 24 chromosomes
def load_genome (genome_filename):
	# read in the genome.fa file using SeqIO
	genome=list(SeqIO.parse(genome_filename, 'fasta'))
	# only take chr1~22 + chrX + chrY
	chrs=[]
	for i in range(22):
		chrs.append('chr'+str(i+1))
	chrs.append('chrX')
	chrs.append('chrY')
	ids=[]
	for i in range(len(genome)):
		ids.append(genome[i].id)
	# put the 24 chrs into chrs_clean
	chrs_clean=[]
	for i in range(len(chrs)):
		chr=chrs[i]
		id_index=ids.index(chr)
		chrs_clean.append(genome[id_index])
	return chrs_clean

# extract sequence based on positions, "strand" specifies which column contains the strandness(+/-)
def seq_extract (chrs_clean, query_filename, strand, length):
	with open(query_filename) as fp:
		out_line=[]
		allids=[x.id for x in chrs_clean]
		for i, line in enumerate(fp):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])
			end=int(line.split('\t')[2])
			if chr in allids:
				index=[m for m in range(24) if chrs_clean[m].id==chr][0]
				sequence=chrs_clean[index].seq[start:end]
				strand=line.split('\t')[strand-1]
				if strand=='-':
					sequence=sequence.reverse_complement()
				if abs(end - start)>=length:
					out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+sequence)
	fp.close()
	f=open(query_filename+'_seqs.fa','w+')
	f.write('\n'.join(str(a) for a in out_line))
	f.close()


# extract sequence based on positions, strandness (+/-) is not specified.
def seq_extract_nostrand (chrs_clean, query_filename, length):
	with open(query_filename) as fp:
		out_line=[]
		allids=[x.id for x in chrs_clean]
		for i, line in enumerate(fp):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])
			end=int(line.split('\t')[2])
			if start < end:
				if chr in allids:
					index=[m for m in range(24) if chrs_clean[m].id==chr][0]
					sequence=chrs_clean[index].seq[start:end]
					if abs(end - start)>=length:
						out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+sequence)
			if start > end:
				if chr in allids:
					index=[m for m in range(24) if chrs_clean[m].id==chr][0]
					sequence=chrs_clean[index].seq[end:start]
					sequence=sequence.reverse_complement()
					if abs(end - start)>=length:
						out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+sequence)
	fp.close()
	f=open(query_filename+'_seqs.fa','w+')
	f.write('\n'.join(str(a) for a in out_line))
	f.close()


# extract L1 5'UTR sequence: before ATGGGG
# 4th column is the strandness
def seq_extract_l1_5utr (chrs_clean, query_filename, length):
	with open(query_filename) as fp:
		out_line=[]
		for i, line in enumerate(fp):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])
			end=int(line.split('\t')[2])
			index=[m for m in range(24) if chrs_clean[m].id==chr][0]
			sequence=chrs_clean[index].seq[start:end]
			strand=str(line.strip().split('\t')[3])
			#print strand
			if strand == '-':
				sequence=sequence.reverse_complement()
			if (end - start)>=length:
				sequence=sequence.upper()
				pos=sequence.find('ATGGGG')
				sequence=sequence[0:pos]
				out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+str(sequence))
	fp.close()
	#print len(out_line)
	f=open(query_filename+'_'+str(length)+'_5UTR'+'_seqs.fa','w+')
	f.write('\n'.join(str(a) for a in out_line))
	f.close()

# this function splits a fasta file into multiple files, each containing one sequence
def split_fasta (fasta_filename):
	fasta=list(SeqIO.parse(fasta_filename, 'fasta'))
	for i in range(len(fasta)):
		sequence=fasta[i].seq
		id=str('>'+fasta[i].id)
		f=open(fasta[i].id+'.fa','w+')
		f.write(id+'\n'+str(sequence))
		f.close()



######################################################################################
######################################################################################
######################################################################################
# main function

def main():
	logging.info('Parsing command line.')
	usage = '%prog [options]'
	parser = optparse.OptionParser(usage=usage)
	parser.add_option(
		'-g','--genome', type='string',
		action='store', dest='genome_filename',
		default='genome/hg38',
		help='Path to the human genome fasta file. Default is \'genome/hg38.fa\' and is most likely not working for you. Please specify.')
	parser.add_option(
		'-q','--query_file', type='string',
		action='store', dest='query_filename',
		default='input',
		help='Path to the query file (in bed format). Default = \'input\' and is most likely not working for you. Please specify.')
	parser.add_option(
		'-s','--strand', type='int',
		action='store', dest='strand',
		help='Column which contains the strandness(+/-) information. Default is None (not considering the strandness).')
	parser.add_option(
		'-l','--length', type='int',
		action='store', dest='length',
		default=0,
		help='Minimum length of sequence to output. Default is 0.')
	(options, args) = parser.parse_args()

	logging.basicConfig(level=logging.INFO,
		format='%(asctime)s %(levelname)s %(message)s')
	#inputs = args[:-1]
	logging.info('Processing input.')
	if options.strand:
		chrs_clean=load_genome(options.genome_filename)
		output_seqs=seq_extract(chrs_clean, options.query_filename, options.strand, options.length)
	else:
		chrs_clean=load_genome(options.genome_filename)
		output_seqs=seq_extract_nostrand(chrs_clean, options.query_filename, options.length)
	logging.info('Final alignment file: {0}'.format(output_seqs))


##############################################
if __name__ == "__main__":
	sys.exit(main())







