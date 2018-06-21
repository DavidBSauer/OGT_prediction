import logging
from matplotlib_venn import venn3
from matplotlib_venn import venn2 
import matplotlib.pyplot as plt

def calc(species_features,species_OGT,species_level_dict):
	shared_species = len([species for species in species_features.keys() if species in species_OGT.keys()]) #find those species with both features and OGTs
	v2 = venn2(subsets = (len(species_features.keys())-shared_species,len(species_OGT.keys())-shared_species,shared_species),set_labels = ('Species w/\nFeature Data', 'Species w/\nOGT Data')) #plot venn diagram of species-(any)features vs species-OGTs
	plt.savefig('./figures/species_feature_OGT.png')
	plt.cla()
	plt.clf()
	plt.close()

	#calculate venn diagram values
	species_list = []
	genomes = []
	tRNAs = []
	rRNAs = []
	ORFs = []
	for species in species_features:
		species_list.append(species)
		for feature in species_features[species]:
			if 'genomic' in feature:
				genomes.append(species)
			elif 'tRNA' in feature:
				tRNAs.append(species)
			elif 'rRNA' in feature:
				rRNAs.append(species)
			elif 'ORF' in feature:
				ORFs.append(species)

	genomes = list(set(genomes))
	tRNAs = list(set(tRNAs))
	rRNAs = list(set(rRNAs))
	ORFs = list(set(ORFs))

	logging.info('the number of species with genomes ='+str(len(genomes)))
	logging.info('the number of species with tRNAs ='+str(len(tRNAs)))
	logging.info('the number of species with rRNAs ='+str(len(rRNAs)))
	logging.info('the number of species with ORFs ='+str(len(ORFs)))

	tRNA_only = len([x for x in tRNAs if not(x in rRNAs) and not(x in ORFs)])
	rRNA_only = len([x for x in rRNAs if not(x in tRNAs) and not(x in ORFs)])
	ORFs_only = len([x for x in ORFs if not(x in tRNAs) and not(x in ORFs)])
	logging.info('the number of species with tRNAs only ='+str(tRNA_only))
	logging.info('the number of species with rRNAs only ='+str(rRNA_only))
	logging.info('the number of species with ORFs only ='+str(ORFs_only))

	tRNA_rRNA = len([x for x in tRNAs if x in rRNAs and not(x in ORFs)])
	tRNA_ORFs = len([x for x in tRNAs if x in ORFs and not(x in rRNAs)])
	rRNA_ORFs = len([x for x in rRNAs if x in ORFs and not(x in tRNAs)])
	logging.info('the number of species with tRNAs and rRNAs but not ORFs ='+str(tRNA_rRNA))
	logging.info('the number of species with tRNAs and ORFs but not rRNA ='+str(tRNA_ORFs))
	logging.info('the number of species with rRNAs and ORFs but not tRNA ='+str(rRNA_ORFs))

	tRNA_rRNA_ORFs = len([x for x in tRNAs if x in rRNAs and x in ORFs])
	logging.info('the number of species with tRNA, rRNAs, and ORFs ='+str(tRNA_rRNA_ORFs))

	#create venn diagram of all species
	v3 =venn3(subsets = (tRNA_only,rRNA_only,tRNA_rRNA,ORFs_only,tRNA_ORFs,rRNA_ORFs,tRNA_rRNA_ORFs),set_labels = ('tRNA only','rRNA only','tRNA rRNA','ORFs only','tRNA ORFs','rRNA ORFs','tRNA rRNA ORFs'))
	plt.savefig('./figures/feature_assignment_all.png')
	plt.cla()
	plt.clf()
	plt.close()

	f = open('./files/species_list_all.txt','w')
	f.write('\n'.join(species_list))
	f.close()

	
	for clade in ['Bacteria','Archaea']:
		species_list = []
		genomes = []
		tRNAs = []
		rRNAs = []
		ORFs = []
		for species in species_features:
			if species in species_level_dict['superkingdom'][clade]:
				species_list.append(species)
				for feature in species_features[species]:
					if 'genomic' in feature:
						genomes.append(species)
					elif 'tRNA' in feature:
						tRNAs.append(species)
					elif 'rRNA' in feature:
						rRNAs.append(species)
					elif 'ORF' in feature:
						ORFs.append(species)

		genomes = list(set(genomes))
		tRNAs = list(set(tRNAs))
		rRNAs = list(set(rRNAs))
		ORFs = list(set(ORFs))
				
		logging.info('the number of '+clade+' species with genomes ='+str(len(genomes)))
		logging.info('the number of '+clade+' species with tRNAs ='+str(len(tRNAs)))
		logging.info('the number of '+clade+' species with rRNAs ='+str(len(rRNAs)))
		logging.info('the number of '+clade+' species with ORFs ='+str(len(ORFs)))


		tRNA_only = len([x for x in tRNAs if not(x in rRNAs) and not(x in ORFs)])
		rRNA_only = len([x for x in rRNAs if not(x in tRNAs) and not(x in ORFs)])
		ORFs_only = len([x for x in ORFs if not(x in tRNAs) and not(x in ORFs)])
		logging.info('the number of '+clade+' species with tRNAs only ='+str(tRNA_only))
		logging.info('the number of '+clade+' species with rRNAs only ='+str(rRNA_only))
		logging.info('the number of '+clade+' species with ORFs only ='+str(ORFs_only))

		tRNA_rRNA = len([x for x in tRNAs if x in rRNAs and not(x in ORFs)])
		tRNA_ORFs = len([x for x in tRNAs if x in ORFs and not(x in rRNAs)])
		rRNA_ORFs = len([x for x in rRNAs if x in ORFs and not(x in tRNAs)])
		logging.info('the number of '+clade+' species with tRNAs and rRNAs but not ORFs ='+str(tRNA_rRNA))
		logging.info('the number of '+clade+' species with tRNAs and ORFs but not rRNA ='+str(tRNA_ORFs))
		logging.info('the number of '+clade+' species with rRNAs and ORFs but not tRNA ='+str(rRNA_ORFs))

		tRNA_rRNA_ORFs = len([x for x in tRNAs if x in rRNAs and x in ORFs])
		logging.info('the number of '+clade+' species with tRNA, rRNAs, and ORFs ='+str(tRNA_rRNA_ORFs))

		v3 =venn3(subsets = (tRNA_only,rRNA_only,tRNA_rRNA,ORFs_only,tRNA_ORFs,rRNA_ORFs,tRNA_rRNA_ORFs),set_labels = ('tRNA only','rRNA only','tRNA rRNA','ORFs only','tRNA ORFs','rRNA ORFs','tRNA rRNA ORFs'))
		plt.savefig('./figures/feature_assignment_'+clade+'.png')
		plt.cla()
		plt.clf()
		plt.close()
		f = open('./files/species_list_'+clade+'.txt','w')
		f.write('\n'.join(species_list))
		f.close()
