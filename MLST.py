import argparse
import os
import sys
import subprocess
import re
from Bio import SeqIO

# use argparse to run through the command line options given
args_parser = argparse.ArgumentParser(description="Program for creating phylogenes the concatenated amino acid sequences of marker genes", epilog="Center for Microbial Oceanography Research and Education")
args_parser.add_argument('-i', '--input_folder', required=True, help='Provide a folder of input sequences.')
args_parser.add_argument('-o', '--output_folder', required=True, help='Provide a folder for output sequences.')
args_parser.add_argument('-d', '--database', required=True, help='Provide a file of the hmm sequences to be matched.')
args_parser = args_parser.parse_args()

# set up object names for input/output/database folders
output_dir = args_parser.output_folder
input_dir = args_parser.input_folder
db = args_parser.database

# Loop through input files given and start HMMER searches
dirs = os.listdir(input_dir)
#print dirs
print "\nFinding proteins encoded in genomes using HMMER..."
for filenames in dirs:
	if filenames.endswith(".faa"):
		cmd = "hmmsearch -E 1e-5 --cpu 3 --tblout " + output_dir+"/"+filenames+".hmmsearch" + " " + db + " " + input_dir+"/"+filenames
		#print cmd,"\n"
		subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
		#os.system(cmd) [Probably better to use subprocess rathern than the typical os.system call

# Loop through and parse HMMER output

dirs = os.listdir(output_dir)
hit_list = {}
score_list = {}
prot_list = {}

for filenames in dirs:
	if filenames.endswith(".hmmsearch"):
		#print output_dir+"/"+filenames
		f = open(output_dir+"/"+filenames, 'r')
		o = open(output_dir+"/"+filenames+".parsed", 'w')
		hit_list = {}
		score_list = {}
		for line in f.readlines():
			if line.startswith("#"):
				list = []
			else:
				#newline = line.replace("\s", "\t")
				newline = re.sub( '\s+', '\t', line)
				list1 = newline.split('\t')
				ids = list1[0]
				hit = list1[2]
				bit_score = list1[5]
				score = float(bit_score)
				#score = "%05d" % (float(bit_score),)
				prot_list[hit] = 1
				if hit in hit_list:
					oldhit_id = hit_list[hit]
					oldhit_score = score_list[hit]
					if score > oldhit_score:
						del hit_list[hit]
						hit_list[hit] = ids

						del score_list[hit]
						score_list[hit] = score	
				else:
					hit_list[hit] = ids
					score_list[hit] = score
		for item in hit_list:
			o.write(item + ' ' + hit_list[item] + ' ' + str(score_list[item]) + '\n')
		o.close()

### Initialize output files for protein sequences
for protid in prot_list:
	#print protid
	o = open(output_dir+'/'+protid+'.hitseqs', 'w')
	o.close()

### Loop through parsed hmmer output and retrieve best hits to the sequences of interest
dirs = os.listdir(output_dir)
mlist = []
for hmm_out in dirs:
	if hmm_out.endswith(".parsed"):
		prot_list2 = {}
		parsed = open(output_dir+'/'+hmm_out, 'r')
		proteins = re.sub('.hmmsearch.parsed', '', hmm_out)
		
		for line in parsed.readlines():
			list1 = line.split(' ')
			prot_id = list1[0]
			seq_id = list1[1]
			seq_score = float(list1[2])

			if prot_id in mlist:
				pass
			else:
				mlist.append(prot_id) 

			if seq_score > 300:
				prot_list2[prot_id] = seq_id
			
		for prot_id in prot_list2:
			for seq_record in SeqIO.parse(input_dir+'/'+proteins, "fasta"):
				if(seq_record.id == prot_list2[prot_id]):
					o = open(output_dir+'/'+prot_id+'.hitseqs', 'a')
					#print seq_record.id
					seq = str(seq_record.seq)
					o.write('>'+seq_record.id+'\n'+seq+'\n')
					o.close()
#print mlist
### Align sequences in the output files using MAFFT
dirs = os.listdir(output_dir)
#print dirs
print "Aligning sequences using MAFFT..."
for filenames in dirs:
	if filenames.endswith(".hitseqs"):
		cmd = "mafft --auto "+output_dir+'/'+filenames+' > '+output_dir+'/'+filenames+'.mafft'
		#print cmd,"\n"
		#subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'w'))
		subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
		#p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
		#out, err = p.communicate()
		#subprocess.call(cmd, shell=True, stdout=f)

### Initialize output files for protein sequences
#for protid in prot_list:
#	#print protid
#	o = open(output_dir+'/'+protid+'.aligned.hitseqs', 'w')
#	o.close()

### Loop through parsed hmmer output and retrieve best hits to the sequences of interest
dirs = os.listdir(output_dir)
o = open(output_dir+'/final.align', 'w')
for hmm_out in dirs:
	if hmm_out.endswith(".parsed"):
		prot_list2 = {}
		parsed = open(output_dir+'/'+hmm_out, 'r')
		#proteins = re.sub('.hmmsearch.parsed', '', hmm_out)
		filename = re.sub('.hmmsearch.parsed', '', hmm_out)
		o.write('>'+filename+'\n')

		for line in parsed.readlines():
			list1 = line.split(' ')
			prot_id = list1[0]
			seq_id = list1[1]
			seq_score = float(list1[2])
			if seq_score > 300:
				prot_list2[prot_id] = seq_id
			#else:
				#prot_list2[prot_id] = 'missing'

		for item in mlist:
			seqs = []
			#o = open(output_dir+'/'+item+'.aligned.hitseqs', 'a')
			
			if item in prot_list2:
				for seq_record in SeqIO.parse(output_dir+'/'+item+'.hitseqs.mafft', "fasta"):
					if seq_record.id == prot_list2[item]:
					#protein_id = prot_list2[item]
						seq = str(seq_record.seq)
						o.write(seq)
						break
			else:
				#o.write('XXXXXXXXXX')
				for seq_record in SeqIO.parse(output_dir+'/'+item+'.hitseqs.mafft', "fasta"):
					length = len(str(seq_record.seq))
					dummy = 'X' * length
					o.write(dummy)
					print filename, item, length
					break
				#for seq_record in SeqIO.parse(output_dir+'/'+item+'.hitseqs.mafft', "fasta"):
				#	representative = seqs[0]
				#	length = len(representative)
				#	dummy = 'X' * length
				#	o.write(dummy+'\n')
			#o.close()
		o.write('\n')
o.close()

print 'Creating phylogeny based on concatenated alignment using FastTree...'
cmd = 'bin/./FastTree '+ output_dir + '/final.align > final.nwk'
subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
print 'Done!\n'










