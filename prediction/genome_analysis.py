from sys import argv
import logging
import csv
import multiprocessing as mp
import os
import genomic
import tRNA
import protein
import ORFs
import external_tools
import rRNA
import operator
import random
from tqdm import tqdm

species_clade = {}
def setup(sc_input):
	global species_clade
	species_clade = sc_input

#Analyze features per genome
def features_per_genome((genome_file,species)):
	result = {}
	result['species']=species
	#create folders and decompress
	genome_file = external_tools.setup((genome_file,species))	
	
	#calculate genomic features
	logging.info('Working on '+genome_file)
	result['genomic'] = genomic.analysis((genome_file,species))
	
	#calculate tRNA features
	(tRNA_test,tRNA_seqs) = external_tools.tRNA((genome_file,species))
	if tRNA_test:
		#if tRNAs were predicted, calculate features
		tRNA_data = tRNA.analysis(tRNA_seqs)
		logging.info('For '+genome_file+' number of tRNAs: '+str(tRNA.number(tRNA_seqs)))
		result['tRNA']=tRNA_data
	else:
		logging.info('Cound not retrieve tRNAs for: '+genome_file)

	#calculate rRNA features
	logging.info('Finding rRNAs for: '+genome_file)
	#run barrnap using both archaea and bacteria hmm models
	domain_results = external_tools.rRNA((genome_file,species))
	#using previous assigned domain
	domain_assigned = species_clade[species]
	if domain_results[domain_assigned]:
		(rRNA_assigned_test,rRNA_assigned_seq) = external_tools.rRNA_seq((genome_file,species),domain_assigned,'assigned')
		if rRNA_assigned_test:		
			rRNA_assigned_data = rRNA.analysis(rRNA_assigned_seq)
			result['rRNA']=rRNA_assigned_data

	#calculate ORF and proteome features
	(ORF_test,ORF_seqs) = external_tools.genemark((genome_file,species))
	if ORF_test:	
		logging.info('Analyzing ORFs for '+genome_file)	
		t_size = result['genomic']['Total Size']
		ORF_data = ORFs.analysis(ORF_seqs,t_size)
		protein_data = protein.analysis(ORF_seqs)
		result['ORF']=ORF_data
		result['protein']=protein_data
	else:
		logging.info('Problem getting ORFs for '+genome_file)	
		
	external_tools.cleanup((genome_file,species))
	return (genome_file,result)

def file_checksum((genome_file,species)):
	return ((genome_file,species),int(os.stat('./genomes/'+species+'/'+genome_file).st_size))



def many_genomes(genomes):
	#sort the genomes by file size to optimize parallelization
	logging.info('sorting genomes by reverse file-size order to improve parallelism')
	print 'analyzing genomes'
	to_analyze = []
	for x in genomes.keys():
		to_analyze.append((x,genomes[x]))	
	p = mp.Pool()
	results = p.map(file_checksum, to_analyze)
	p.close()
	results = {x[0]:x[1] for x in results}
	sortedgenomes = [x[0] for x in sorted(results.items(),key=operator.itemgetter(1))]
	sortedgenomes.reverse() #largest to smallest for analysis, comment out for trouble shooting

	#make folders for output
	if not(os.path.isdir('./output')):	
		os.mkdir('./output')
	if not(os.path.isdir('./output/genomes')):	
		os.mkdir('./output/genomes')

	for species in [x[1] for x in sortedgenomes]:
		#make a folder based on the genome name 
		if not(os.path.isdir('./output/genomes/'+species)):	
			os.mkdir('./output/genomes/'+species)

	
	#calculate features in parallel
	p = mp.Pool()
	results = p.map(features_per_genome, sortedgenomes)
	p.close()
	'''
	#calculate single thread for troubleshooting
	results =[]
	for x in tqdm(sortedgenomes,unit='genome'):
		results.append(features_per_genome(x))
	'''
	results = {x[0]:x[1] for x in results}
	
	remod_results ={genome:{} for genome in results.keys()}
	#write out the results by feature class
	for feature_class in ['genomic','tRNA','rRNA','ORF','protein']:
		all_params = []
		for genome in results.keys():
			if feature_class in results[genome].keys():
				all_params = results[genome][feature_class].keys()
				break

		g= open('Genome_'+feature_class+'_features.txt','w')
		g.write('Genome\tspecies\t'+'\t'.join([feature_class+' '+x for x in all_params])+'\n')
		for genome in results.keys():
			if feature_class in results[genome].keys():
				g.write(genome+'\t'+results[genome]['species']+'\t')
				for param in all_params:
					g.write(str(results[genome][feature_class][param])+'\t')
					remod_results[genome][feature_class+' '+param] = results[genome][feature_class][param]
				g.write('\n')
		g.close()
	for genome in results.keys():
		remod_results[genome]['species'] = results[genome]['species']
	return remod_results
