from sys import argv
import logging
import csv
import multiprocessing as mp
import random
import genomic
import tRNA
import protein
import ORFs
import external_tools
import rRNA
import operator
import os

script, genome_species_file, species_clade_file = argv

logging.basicConfig(filename=str('feature_calculation.log'), level=logging.INFO)

#read in the species-clade file. save as dict species->clade
logging.info("trying to open species-taxon file")
infile = open(species_clade_file,'r')
species_clade = {}
for line in infile.readlines():
	working = line.split('\t')
	if not(working[1] == 'None'):
		species_clade[working[0]] = working[1].strip() 
infile.close()    
logging.info("found "+str(len(species_clade.keys()))+" species with taxonomic information")

#read in the genomes file. save a dict of genome->species
logging.info("trying to open genomes file")
infile = open(genome_species_file,'r')
reader = csv.reader(infile,delimiter='\t')
genomes = dict((str(rows[0]),str(rows[1])) for rows in reader)
infile.close()    
logging.info("found "+str(len(genomes.keys()))+" genomes")

#sort the genomes by file size to optimize parallelization
to_analyze = [(x,genomes[x]) for x in genomes.keys()]
random.shuffle(to_analyze)

#make folders for output
if not(os.path.isdir('./output')):	
	os.mkdir('./output')
if not(os.path.isdir('./output/genomes')):	
	os.mkdir('./output/genomes')

for species in [x[1] for x in to_analyze]:
	#make a folder based on the genome name 
	if not(os.path.isdir('./output/genomes/'+species)):	
		os.mkdir('./output/genomes/'+species)

logging.info('Analyzing Genomic features')
#Analyze features per genome
def features_per_genome(inputs):
	(genome_file,species) = inputs
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
	#attempt to predict domain using 16S rRNA sequences
	(pred_test,domain_pred) = external_tools.classify((genome_file,species))
	#calculate rRNA features based on predicted domain
	if pred_test:
		result['rRNA_domain']=domain_pred	
		(rRNA_pred_test,rRNA_pred_seq) = external_tools.rRNA_seq((genome_file,species),domain_pred,'pred')
		if rRNA_pred_test:		
			rRNA_pred_data = rRNA.analysis(rRNA_pred_seq)
			result['rRNA_pred']=rRNA_pred_data
	#using previous assigned domain
	global species_clade
	domain_assigned = species_clade[species]
	if domain_results[domain_assigned]:
		(rRNA_assigned_test,rRNA_assigned_seq) = external_tools.rRNA_seq((genome_file,species),domain_assigned,'assigned')
		if rRNA_assigned_test:		
			rRNA_assigned_data = rRNA.analysis(rRNA_assigned_seq)
			result['rRNA_assigned']=rRNA_assigned_data

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

#calculate features in parallel
'''
p = mp.Pool()
results = p.map(features_per_genome, to_analyze)
p.close()
'''
results =[]
for x in to_analyze:
	results.append(features_per_genome(x))

results = {x[0]:x[1] for x in results}

#write out the results by feature class
for feature_class in ['genomic','tRNA','rRNA_pred','rRNA_assigned','ORF','protein']:
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
			g.write('\n')
	g.close()

#write out barrnap assignment of domain
h= open('Genome_barrnap_assignment.txt','w')
for genome in results.keys():
	if 'rRNA_domain' in results[genome].keys():
		h.write(genome+'\t'+results[genome]['rRNA_domain']+'\n')
h.close()

logging.info('Exiting normally')
