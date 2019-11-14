
#This Python script requires Biopython 1.51 or later
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Alphabet import generic_dna
import itertools
import gzip
from Bio.Seq import Seq
import pandas as pd
import sys
import numpy as np
import shutil
import os
import tqdm
#import subprocess

def get_len(in_file):
	def blocks(files, size=65536):
		while True:
			b = files.read(size)
			if not b: break
			yield b
	
	with gzip.open(in_file,"rU") as f:
		res = sum(bl.count("\n") for bl in blocks(f))
		res = float(res)/3   
		return res

def find_barcode(barcode = '', in_seq=''):
	barcode = Seq(barcode, generic_dna)
	rc_barcode = barcode.reverse_complement()
	in_seq = Seq(in_seq, generic_dna)
	if barcode in in_seq:
		return 1, 'f'
	elif rc_barcode in in_seq:
		return 1, 'r'
	else:
		return 0, 'n'

def remove_barcode(seq, quol, barcode_orient):
	if barcode_orient == 'f':
		start = seq.find(f_barcode)+len(f_barcode)
		seq = seq[start:]
		quol = quol[start:]
	if barcode_orient == 'r':
		start = seq.find(r_barcode)
		seq = seq[:start]
		quol = quol[:start]
	return   seq,  quol  
         
         
def parse(file_name, f_seq_barcode, r_seq_barcode, inFileLength):
	print('Search for:', f_seq_barcode, r_seq_barcode)
	print('Setup results files')
	file_r1 = file_name+"1.fq.gz"
	file_r2 = file_name+"2.fq.gz"
	print (inFileLength)
	ff_barcode_1 = open(file_name+"1_ff_barcode.fq",'w')
	ff_barcode_2 = open(file_name+"2_ff_barcode.fq",'w')
	
	fr_barcode_1 = open(file_name+"1_fr_barcode.fq",'w')
	fr_barcode_2 = open(file_name+"2_fr_barcode.fq",'w')
	
	
	rr_barcode_1 = open(file_name+"1_rr_barcode.fq",'w')
	rr_barcode_2 = open(file_name+"2_rr_barcode.fq",'w')
	
	rf_barcode_1 = open(file_name+"1_rf_barcode.fq",'w')
	rf_barcode_2 = open(file_name+"2_rf_barcode.fq",'w')
	
	count = 0
	count_ff = 0
	count_fr = 0
	count_rr = 0
	count_rf = 0
	
	r1_iter = FastqGeneralIterator(gzip.open(file_r1,"rt"))
	r2_iter = FastqGeneralIterator(gzip.open(file_r2,"rt"))
	for (r1_id, r1_seq, r1_q), (r2_id, r2_seq, r2_q) in tqdm.tqdm(zip(r1_iter, r2_iter), total=int(inFileLength), miniters=500000):
	#for (r1_id, r1_seq, r1_q), (r2_id, r2_seq, r2_q) in zip(r1_iter, r2_iter):
		#test we are looking at paired mate
		#print( f_id,r_id)
		#assert r1_id.split(' ')[0] == r2_id.split(' ')[0]
		count += 1

		#look at forward barcode
		found_r1, orient_r1 =   find_barcode(barcode = f_seq_barcode, in_seq=r1_seq)
		found_r2, orient_r2 =   find_barcode(barcode = f_seq_barcode, in_seq=r2_seq)

		if found_r1 == 1 and orient_r1 == 'f':
			ff_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			ff_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_ff+=1
		if found_r2 == 1 and orient_r2 == 'f':
			ff_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			ff_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_ff+=1
		
		if found_r1 == 1 and orient_r1 == 'r':
			fr_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			fr_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_fr+=1
		if found_r2 == 1 and orient_r2 == 'r':
			fr_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			fr_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_ff+=1

		#look at reverse barcode
		found_r1, orient_r1 =   find_barcode(barcode = r_seq_barcode, in_seq=r1_seq)
		found_r2, orient_r2 =   find_barcode(barcode = r_seq_barcode, in_seq=r2_seq)

		if found_r1 == 1 and orient_r1 == 'f':
			rf_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			rf_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_rf+=1
		if found_r2 == 1 and orient_r2 == 'f':
			rf_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			rf_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_rf+=1
		
		if found_r1 == 1 and orient_r1 == 'r':
			rr_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			rr_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_rr+=1
		if found_r2 == 1 and orient_r2 == 'r':
			rr_barcode_1.write("@%s\n%s\n+\n%s\n" % (r1_id, r1_seq, r1_q))
			rr_barcode_2.write("@%s\n%s\n+\n%s\n" % (r2_id, r2_seq, r2_q))
			count_rr+=1

		
		
	ff_barcode_1.close()
	ff_barcode_2.close()
	fr_barcode_1.close()
	fr_barcode_2.close()
	
	
	rr_barcode_1.close()
	rr_barcode_2.close()
	rf_barcode_1.close()
	rf_barcode_2.close()
	
	def gzip_fastq(infile):
		with open(infile, 'rb') as f_in:
			with gzip.open(infile+'.gz', 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
				os.remove(infile)
	
	gzip_fastq(ff_barcode_1.name)
	gzip_fastq(ff_barcode_2.name)
	
	gzip_fastq(fr_barcode_1.name)
	gzip_fastq(fr_barcode_2.name)
	
	gzip_fastq(rr_barcode_1.name)
	gzip_fastq(rr_barcode_2.name)
	
	gzip_fastq(rf_barcode_1.name)
	gzip_fastq(rf_barcode_2.name)

	print('all:',count, 'ff:',count_ff, 'fr:', count_fr, 'rf:', count_rf, 'rr:', count_rr)
	
if __name__ == '__main__':

	#barcode_dictionary ={
	#'CosLib':{'f':'CTCTTAAAAGCATCATGTCT', 'r':'ACTAGTTCTAGAGCGGCCGC'},
	#'OElib':{'f':'GATAGAGTGGTACCGGCCGG', 'r':'CAATGATAGAGTGGCCGGCC'}
	#}
	barcode_dictionary ={
	'CosLib':{'f':'CTCTTAAAAGCATCATGTCT', 'r':'ACTAGTTCTAGAGCGGCCGC'},
	'OElib':{'f':'GATAGAGTGGTACCGGCCGG', 'r':'CAATGATAGAGTGGCCGGCC'}
	}
	r_barcode = barcode_dictionary[sys.argv[1]]['f']
	f_barcode = barcode_dictionary[sys.argv[1]]['r']
	parse(sys.argv[2], f_barcode, r_barcode, sys.argv[3])
