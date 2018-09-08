import logging
import time
from Bio import Entrez
from tqdm import tqdm

import datetime
from sys import argv
logging.basicConfig(filename=str('clade_retriever.log'), level=logging.INFO)
logging.info('accessing NCBI-Taxonomy on: '+datetime.datetime.utcnow().strftime('%m/%d/%Y'))

species_list_input = argv[1]
Entrez.email = argv[2]
logging.info('species list: '+species_list_input)
logging.info('email address: '+argv[2])

def clade(species_list):
	taxa_info = {}
	for name in tqdm(species_list,unit='species'):
		if len(name.split('_'))>= 2:
			species=name.lower().replace('_','+')
			try: 
				search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
			except:
				logging.info('Trouble getting the taxid for: '+species)
			else:
				record = Entrez.read(search)
				if len(record['IdList'])>0:
					taxid = record['IdList'][0]
					try:
						search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
					except:
						logging.info('Trouble getting the taxonomy info for: '+taxid)
					else:
						result = Entrez.read(search)
						if len(result)>0:
							lineage = result[0]['LineageEx']
							taxa_info[name] = {}
							for level in result[0]['LineageEx']:
								if level['Rank'] == 'superkingdom':
									taxa_info[name]['superkingdom'] = level['ScientificName']
								if level['Rank'] == 'class':
									taxa_info[name]['class'] = level['ScientificName']
								if level['Rank'] == 'phylum':
									taxa_info[name]['phylum'] = level['ScientificName']
								if level['Rank'] == 'order':
									taxa_info[name]['order'] = level['ScientificName']
								if level['Rank'] == 'family':
									taxa_info[name]['family'] = level['ScientificName']

	g = open('species_taxonomic.txt','w')
	g.write('\t'.join(['species','superkingdom','phylum','class','order','family'])+'\n')
	for species in taxa_info.keys():
		g.write(species)
		for level in ['superkingdom','phylum','class','order','family']:
			if level in taxa_info[species].keys():
				g.write('\t'+taxa_info[species][level])
			else:
				g.write('\tNone')
		g.write('\n')
	g.close()
	return taxa_info

f= open(species_list_input,'r')
lines = [line for line in f.readlines() if len(line.split())>0]
f.close()
species_list_input = [line.split()[0].strip() for line in lines]
logging.info('the number of input species: '+str(len(species_list_input)))
logging.info('the number of species with retrieved clade info: '+str(len(clade(species_list_input).keys())))

