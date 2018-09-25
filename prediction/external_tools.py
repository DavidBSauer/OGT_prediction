import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_rna, generic_dna
from Bio.Seq import Seq
import shutil
import numpy as np
import subprocess
import logging
logger = logging.getLogger('feature_calculation')
from BCBio import GFF
import gzip
import tarfile

#load external tool commands
commands = {}
f = open('external_tools.txt','r')
for line in f.readlines():
	commands[line.split()[0].strip()]=line.split()[1].strip()
f.close()

#record version
def versions():
	global commands
	p = subprocess.Popen([commands['tRNAscan-SE']+' -h'],shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	logger.info('tRNAscan-SE version info: '+err.decode('utf-8').split('\n')[1].strip())
	p = subprocess.Popen([commands['bedtools']+' --version'],shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	logger.info('bedtools version info: '+out.decode('utf-8').strip())
	p = subprocess.Popen([commands['barrnap']+' --version'],shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	logger.info('barrnap version info: '+err.decode('utf-8').strip())
	p = subprocess.Popen([commands['genemark']+' --version'],shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	logger.info('genemark version info: '+out.decode('utf-8').strip())
	logger.info('Numpy version: '+np.__version__)
	import Bio
	logger.info('Biopython version: '+Bio.__version__)
	del(Bio)
	import sys
	logger.info('Python version: '+sys.version)
	del(sys)
	import platform
	logger.info('Platform: '+platform.platform())
	del(platform)

versions()

def make_tarfile(inputs):
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	with tarfile.open('./output/genomes/'+species+'/'+folder+'.tar.gz', "w:gz") as tar:
		tar.add('./output/genomes/'+species+'/'+folder, arcname=os.path.basename('./output/genomes/'+species+'/'+folder))

#decompress genome
def setup(inputs):
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-2])
	os.mkdir('./output/genomes/'+species+'/'+folder)
	#uncompress genome
	input_file = gzip.open('./genomes/'+species+'/'+genome_file,'r')
	g = open('./output/genomes/'+species+'/'+folder+'/'+'.'.join(genome_file.split('.')[:-1]),'w')
	for line in input_file:	
		g.write(line.decode('ascii').upper())
	g.close()
	input_file.close()
	return '.'.join(genome_file.split('.')[:-1])

#clean up by decompress genome and compressing genemark and barrnap files
def cleanup(inputs):
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	#remove the uncompressed genome file
	os.remove('./output/genomes/'+species+'/'+folder+'/'+genome_file)
	#compress the results
	make_tarfile((genome_file,species))
	#remove uncompressed results
	shutil.rmtree('./output/genomes/'+species+'/'+folder)

#find tRNAs with tRNAScan-SE
def tRNA(inputs):
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	#run tRNAscan
	#settings: -G General, -q Quiet, -b brief output, -o output_file
	global commands
	command = commands['tRNAscan-SE']+' -G -q --brief --output ./output/genomes/'+species+'/'+folder+'/trnascan_gff.txt --struct ./output/genomes/'+species+'/'+folder+'/trnascan_str.txt ./output/genomes/'+species+'/'+folder+'/'+genome_file
	p = subprocess.Popen([command],shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	out = out.decode('utf-8')
	err = err.decode('utf-8')
	if err == '': #catching errors for genomes that fail in genemark analysis
		tRNAs = []
		#only proceed if the tRNAscan-SE results file exists
		if os.path.isfile('./output/genomes/'+species+'/'+folder+'/trnascan_str.txt'):
			f = open('./output/genomes/'+species+'/'+folder+'/trnascan_str.txt','r')
			seqs =[]
			for line in f.readlines():
				if line.startswith('Seq:'):
					new_seq = line.split()[1] 				
					seqs.append(SeqRecord(Seq(str(new_seq),generic_rna),'tRNA_'+str(len(seqs)),'',''))
			f.close()
			SeqIO.write(seqs,'./output/genomes/'+species+'/'+folder+'/trnascan_result.fa','fasta')
			return (True,SeqIO.index('./output/genomes/'+species+'/'+folder+'/trnascan_result.fa','fasta'))
		else:
			#the tRNA file did not exist
			return (False,None)
	else:
		#tRNAscan-SE produced errors
		logger.info('error with tRNAscan for '+genome_file+' with a message of \n'+err)
		return (False,None)

#using GeneMark ORFfinder
def genemark(inputs):
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	os.mkdir('./output/genomes/'+species+'/'+folder+'/genemark')
	os.chdir('./output/genomes/'+species+'/'+folder+'/genemark/')
	#run gmsn.pl
	#settings: --prok prokaryotic, --combine model parameters
	global commands
	command = commands['genemark']+' --prok --combine --gm --fnn ../'+genome_file
	p = subprocess.Popen([command],shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	out = out.decode('utf-8')
	err = err.decode('utf-8')	
	if err == '': #catching errors for genomes that fail in genemark analysis
		#return True if no errors produced by gmhmmp
		os.chdir('../../../../../')
		if os.path.isfile('./output/genomes/'+species+'/'+folder+'/genemark/'+genome_file+'.fnn'):
			return (True,SeqIO.index('./output/genomes/'+species+'/'+folder+'/genemark/'+genome_file+'.fnn','fasta'))
		else:
			logger.info('could not find predicted ORFs for '+genome_file)
			return (False,None)
	else:
		#gmsn.pl produced and error		
		logger.info('error on '+genome_file+' gmsn.pl step with a message of\n'+err)
		os.chdir('../../../../../')
		return (False,None)

#predict rRNA sequences using barrnap
def rRNA(inputs):
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	os.mkdir('./output/genomes/'+species+'/'+folder+'/barrnap')
	domain_results = {}
	#calculate rRNAs using bacterial hmm
	#settings: --quiet Quiet, --threads number_of_threads, --kingdom kingdom_hmm
	global commands
	command = commands['barrnap']+' --quiet --threads 1 --kingdom bac ./output/genomes/'+species+'/'+folder+'/'+genome_file+' > ./output/genomes/'+species+'/'+folder+'/barrnap/Bacteria.txt'
	p = subprocess.Popen(command,shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	out = out.decode('utf-8')
	err = err.decode('utf-8')
	if err == '':
		domain_results['Bacteria'] = True
	else:
		logger.info('error with barrnap bacterial for '+genome_file+' with a message of\n'+err)
		domain_results['Bacteria'] = False
	#calculate rRNAs using archaea hmm
	command = commands['barrnap']+' --quiet --threads 1 --kingdom arc ./output/genomes/'+species+'/'+folder+'/'+genome_file+' > ./output/genomes/'+species+'/'+folder+'/barrnap/Archaea.txt'
	p = subprocess.Popen(command,shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	out = out.decode('utf-8')
	err = err.decode('utf-8')
	if err == '':
		domain_results['Archaea'] = True
	else:
		logger.info('error with barrnap archaea for '+genome_file+' with a message of\n'+err)
		domain_results['Archaea'] = False
	return domain_results
	
#classify based on rRNA sequences
def classify(inputs):
	(genome_file,species) = inputs
	#considering only 16S rRNA sequences. Test both archaea and eukaryotic HMM, take best hit
	#if a tie, this would presume bacterial
	folder = '_'.join(genome_file.split('.')[:-1])
	clade_scores = {}
	ref_recs = SeqIO.to_dict(SeqIO.parse('./output/genomes/'+species+'/'+folder+'/'+genome_file,'fasta'))
	for x in ['Archaea','Bacteria']:
		scores = []
		for b in GFF.parse('./output/genomes/'+species+'/'+folder+'/barrnap/'+x+'.txt',target_lines=1000, base_dict=ref_recs):
			for feature in b.features:
				if feature.qualifiers['product'][0]=='16S ribosomal RNA':
					scores.append(float(feature.qualifiers['score'][0]))
		if len(scores)>0:
			clade_scores[np.mean(scores)]=x
	if len(clade_scores.keys())>0:
		#found 16S rRNA sequences, return best scoring hmm as barrnap assigned domain
		return (True,clade_scores[min(clade_scores.keys())])
	else:
		#could not find any 16S rRNA sequences for classification
		return (False,None)


#calculate rRNA sequences 
def rRNA_seq(inputs):
	(genome_file,species,domain,method) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	#read in the barrnap results file based on the provide domain assignment
	f = open('./output/genomes/'+species+'/'+folder+'/barrnap/'+domain+'.txt','r')
	g = open('./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_'+method+'_16S.txt','w')
	for line in f.readlines():
		if 'Name=16S_rRNA;product=16S ribosomal RNA' in line.split('\t')[-1].strip():
			g.write(line)
	f.close()
	g.close()
	global commands
	command = commands['bedtools']+' getfasta -fi ./output/genomes/'+species+'/'+folder+'/'+genome_file+' -bed ./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_'+method+'_16S.txt -s -name -fo ./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_results_'+method+'.fa'
	p = subprocess.Popen(command,shell=True,executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	out,err = p.communicate()
	out = out.decode('utf-8')
	err = err.decode('utf-8')
	if os.path.isfile('./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_results_'+method+'.fa'):
		#handle if there are multiple rRNA predicted
		f = open('./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_results_'+method+'.fa','r')
		g = open('./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_results_'+method+'_renamed.fa','w')
		count =1
		for line in f.readlines():
			if line[0] == '>':
				g.write(line.replace('>','>'+str(count)+'_'))
				count = count+1
			else:
				g.write(line)
		f.close()
		g.close()
		
		seqs = SeqIO.index('./output/genomes/'+species+'/'+folder+'/barrnap/barrnap_results_'+method+'_renamed.fa','fasta')
		if len(seqs) >0:
			return (True,seqs)
		else:
			#did not predict any rRNA	
			logger.info('did not predict any rRNA for '+genome_file)
			return (False,None)
	else:
		#bedtools produced an error
		logger.info('error in bedtools for '+genome_file+' with a message of\n'+err)
		return (False,None)

