import urllib
from urllib.request import urlopen
import time
import sys
import logging
import os
from tqdm import tqdm
import requests
logging.basicConfig(filename=str('genome_retriever.log'), level=logging.INFO)

ref_file, setting, release = sys.argv[1:4]

logging.info('Reference file: '+ref_file)
logging.info('Setting: '+setting)
logging.info('Using Ensembl release: '+release) #the paper used release 40 but this allow for accessing updated releases

if not(setting in ['IN','NOT_IN']):
	print("Setting must be 'IN' or 'NOT_IN'. Quitting")
	logging.info('Improper setting. Quitting')
	sys.exit()
if setting == 'IN':
	logging.info('Retrieving genomes of species in the provided file.')
else:
	logging.info('Retrieving genomes of species NOT in the provided file.')

logging.info("trying to open species-taxon file")
infile = open(ref_file,'r')
species_list = []
infile.readline()
for line in infile.readlines():
	working = line.split()
	if ((len(working)>0) and (not(working[0].strip() == ''))):
		species_list.append(working[0].strip()) 
infile.close()    
logging.info("found "+str(len(species_list))+" species")
properly_formed = [species for species in species_list if (not(species[0]=='_') and not(species[-1]=='_') and not(species.split('_')[-1] == 'sp.') and len(species.split('_'))==2)]
logging.info('found '+str(len(properly_formed))+' properly formed species names')

print('finding all valid genomes')
  
addresses = {}
root ="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-"+release
retrieve_all= root+"/species_EnsemblBacteria.txt"
try:
	response = urllib.request.urlopen(retrieve_all)
	data = response.read()
	decoded = data.decode('utf-8').split('\n')
except:
	logging.info('could not open the Ensembl bacteria ftp list')
else:
	first_line= decoded[0].split('\t')
	for line in decoded[1:]:
		if not(line.strip() == ''):
			working =line.split('\t')
			x={first_line[pos]:working[pos] for pos in range(0,len(first_line),1)}
			if len(x['species'].split('_'))>=2:
				species ='_'.join(x['species'].split('_')[0:2]).lower()
				(genus,species_name) = species.split('_')
				if not(genus in ['','candidate','uncultured','unidentified','synthetic','candidatus','bacterium','marine']) and not(species_name in ['sp','group','candidatus','bacterium','proteobacterium','endosymbiont','archaeon','cluster','producing','gamma']):
					if not 'bacterium' in species_name:
						if setting == 'IN':
							if species in properly_formed:
								root_address =root+'/fasta/'
								collection = '_'.join(x['core_db'].split('_')[0:3])
								full_name = x['species']
								directory = root_address+collection+'/'+full_name+'/dna/'
								addresses[directory]=species
						else:
							if not(species in properly_formed):
								root_address =root+'/fasta/'
								collection = '_'.join(x['core_db'].split('_')[0:3])
								full_name = x['species']
								directory = root_address+collection+'/'+full_name+'/dna/'
								addresses[directory]=species

logging.info('Number of genomes to retrieve: '+str(len(addresses.keys())))
logging.info('Number of species to retrieve: '+str(len(list(set(addresses.values())))))

if not os.path.exists('./genomes'):
   os.makedirs('./genomes')

def download_file(genome_addr,species):
	try:
		response = urlopen(genome_addr,timeout=60)
		page = response.read().decode('utf-8')
		page = page.split('\n')
		files = [x.strip().split()[-1] for x in page if not(x.strip() == '')]
		g_file = [x for x in files if x.split('.')[-4:] == ['dna','toplevel','fa','gz']]
		if len(g_file)>0:
			g_file = g_file[0]
			if not os.path.exists('./genomes/'+species):
				os.makedirs('./genomes/'+species) 
			with urlopen(genome_addr+g_file, timeout=60) as response, open('./genomes/'+species+'/'+g_file,'wb') as f: 
				the_file = response.read()
				f.write(the_file)
				f.close()
		else:
			return (False,None)
	except:
		return (False,None)
	else:
		return (True,g_file)

print('downloading genomes')
retrieved = {}
t0 = time.time()
for genome_addr in tqdm(addresses.keys(),unit='genome'):
	species = addresses[genome_addr]
	#download the file
	(test,addr) = download_file(genome_addr,species)
	if test:
		retrieved[addr] = species
	else:
		logging.info('problem downloading '+genome_addr)
	#avoid abusing Ensembl server. maximum of one genome per second
	if time.time()-t0 < 1:
		time.sleep(1-(time.time()-t0))
	t0= time.time()

logging.info('of '+str(len(addresses.keys()))+' expected genomes, downloaded '+str(len(retrieved.keys())))
logging.info('from '+str(len(list(set(retrieved.values()))))+' species')

g = open('genomes_retrieved.txt','w')
for genome in retrieved.keys():
	g.write(genome+'\t'+retrieved[genome]+'\n')
g.close()

g = open('species_retrieved.txt','w')
for species in list(set(retrieved.values())):
	g.write(species+'\n')
g.close()
