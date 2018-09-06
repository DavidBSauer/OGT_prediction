from urllib.request import urlopen
import time
from sys import argv
import logging
import csv
import os
import random
from tqdm import tqdm

script, ref_file = argv

if not os.path.exists('./genomes'):
   os.makedirs('./genomes')

 
logging.basicConfig(filename=str('genome_retriever.log'), level=logging.INFO)
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
  
addresses ={}
root_addr= 'ftp://ftp.ensemblgenomes.org/pub/bacteria/release-35/fasta/'
time.sleep(0.1)
try:
	response = urlopen(root_addr, timeout=60)
except:
	#catch if in wrong database
	logging.info('could not open the ftp page')
else:
	the_page = response.read().decode('utf-8')
	the_page = the_page.split('\n')
	collections = [page.split()[-1].strip() for page in the_page if not(page == '')]
	for collection in tqdm(collections,unit='collection'):
		addr = root_addr+collection
		time.sleep(0.1)
		try:
			response = urlopen(addr, timeout=60)
		except:
			#catch if in wrong database
			logging.info('could not open the ftp page for collection: '+collection)
		else:
			the_page = response.read().decode('utf-8')
			the_page = the_page.split('\n')
			genomes = [page.split()[-1].strip() for page in the_page if not(page == '')]
			for genome in genomes:
				if '_'.join(genome.split('.')[0].split('_')[0:2]).lower() in properly_formed:
					genome_addr = addr+'/'+genome+'/dna'
					addresses[genome_addr] = '_'.join(genome.split('.')[0].split('_')[0:2]).lower()

def genome_checksum(genome_addr):
	checksum_file = genome_addr+'/CHECKSUMS'
	try:
		response = urlopen(checksum_file,timeout=60)
	except:
		return (False,None,None)
	else:
		page = response.read().decode('utf-8')
		page = page.split('\n')
		files = {x.split()[-1]:x.split()[0] for x in page if not(x == '')}
		for file in files.keys():
			if len(file.split('.'))>1 and (file.split('.')[-3] == 'toplevel' and file.split('.')[-4] == 'dna'):
				return (True,int(files[file]),file)
		return (False,None,None)

def file_checksum(genome_addr,g_file):
	if g_file != None:
		try:
			species = '_'.join(g_file.split('.')[0].split('_')[0:2]).lower()
			if not os.path.exists('./genomes/'+species):
				os.makedirs('./genomes/'+species) 
			with urlopen(genome_addr+'/'+g_file, timeout=60) as response, open('./genomes/'+species+'/'+g_file,'wb') as f: 
				the_file = response.read()
				f.write(the_file)
				f.close()
		except:
			return (False,None,None)
		else:

			return (True,int(os.popen('sum ./genomes/'+species+'/'+g_file).read().split()[0]),'./genomes/'+species+'/'+g_file)
		return (False,None,None)
	else:
		return (False,None,None)


print('downloading all valid genomes')
retrieved_list ={}
retrieved = 0
logging.info('Number of genomes to retrieve: '+str(len(addresses)))
'''
#retrieve only one genome per species, comment out for final run
old_addresses = addresses
addresses ={}
species =[]
for genome_addr in old_addresses.keys():
	if not(old_addresses[genome_addr] in species):
		addresses[genome_addr]=old_addresses[genome_addr]
'''
for genome_addr in tqdm(addresses.keys(),unit='genome'):
	if random.random() < 1.001: #use <1 (0.001-0.0005) for toubleshooting smaller datasets
		species = addresses[genome_addr]
		(g_test,g_checksum,g_name) = genome_checksum(genome_addr)
		(f_test,f_checksum,f_file) = file_checksum(genome_addr,g_name)
		for z in range(0,10,1):
			if g_checksum != f_checksum or (not(g_test) or not(f_test)):
				print('using checksum loop')
				time.sleep(60)
				(g_test,g_checksum,g_name) = genome_checksum(genome_addr)
				(f_test,f_checksum,f_file) = file_checksum(genome_addr,g_name)
		if g_checksum == f_checksum and g_test:
			retrieved = retrieved+1
			f_file = f_file.split('/')[-1]
			if species in retrieved_list.keys():
				retrieved_list[species].append(f_file)
			else:
				retrieved_list[species] =[f_file]
		else:
			logging.info('Cound not get an checksumed genome for '+genome_addr.split('/')[-2])
			if f_test:
				if os.path.isfile(f_file):
					os.remove(f_file)						

logging.info('of '+str(len(addresses))+' expected genomes, downloaded '+str(retrieved))
logging.info('from '+str(len(retrieved_list.keys()))+' species')

g = open('genomes_retrieved.txt','w')
for x in retrieved_list.keys():
	for y in retrieved_list[x]:
		g.write(y+'\t'+x+'\n')
g.close()
	
