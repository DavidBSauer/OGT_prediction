#take in species-OGT list, species-domain list, and species-feature dictionary
#calculate regression for each feature
#return dict of [domain][feature] = r-value
import logging
logger = logging.getLogger('regression')
import matplotlib
logger.info('Matplotlib version: '+matplotlib.__version__)
matplotlib.use('Agg')
matplotlib.rcParams['axes.linewidth'] = 3
import matplotlib.pyplot as plt
import scipy
logger.info('SciPy version: '+scipy.__version__)
del scipy
from scipy.stats import pearsonr
import numpy as np
logger.info('NumPy version: '+np.__version__)
import os
from matplotlib.ticker import MaxNLocator
import multiprocessing as mp

params = {
   'axes.labelsize': 12,
   'font.size': 12,
   'legend.fontsize': 24,
   'xtick.labelsize': 24,
   'ytick.labelsize': 24,
   'text.usetex': False,
   'figure.figsize': [12, 8]
   }
plt.rcParams.update(params)

if not os.path.exists('./files/individual_features'):
   os.makedirs('./files/individual_features')
if not os.path.exists('./figures/individual_features'):
   os.makedirs('./figures/individual_features')


def r_calc(inputs):
	(feature,species_features,species_OGTs,rank,clade,valid_species,unit) = inputs
	if not os.path.exists('./figures/individual_features/'+rank):
	   os.makedirs('./figures/individual_features/'+rank)
	if not os.path.exists('./files/individual_features/'+rank):
	   os.makedirs('./files/individual_features/'+rank)
	if not os.path.exists('./figures/individual_features/'+rank+'/'+clade):
	   os.makedirs('./figures/individual_features/'+rank+'/'+clade)
	if not os.path.exists('./files/individual_features/'+rank+'/'+clade):
	   os.makedirs('./files/individual_features/'+rank+'/'+clade)

	feature_r = None #catch if a problem in linear regression
	values = []
	plt.cla()
	plt.clf()
	plt.close()
	x_values =[]
	y_values =[]
	#read in all feature values when available
	for species in valid_species:
		if (species in species_features.keys()) and (feature in species_features[species].keys()): #ensure species has feature of interest calculated
			y_values.append(species_features[species][feature])
			x_values.append(species_OGTs[species])
	if (len(x_values)>2) and ((max(y_values)-min(y_values))>0) and ((max(x_values)-min(x_values))>0): #make sure enough values for linear regression (and appropriate values so the regression does not fail)
		#for each domain calculate the correlation between the feature and OGT
		feature_r = pearsonr(x_values,y_values)[0]
	
		#plot the regression for all data
		plt.plot(x_values,y_values,'.',markersize=18,color='#222299',alpha=1)
		#plt.xlim(-10,110)
		#plt.xticks([0,25,50,75,100])
		plt.tick_params(axis='x',which='both',bottom=False,top=False)
		plt.tick_params(axis='y',which='both',left=False,right=False)
		plt.title(feature+'\nr='+str(pearsonr(x_values,y_values)[0])[0:6])
		plt.rc('font')
		plt.ylabel(feature)
		plt.xlabel('Measured '+unit,fontsize=25)
		plt.savefig('./figures/individual_features/'+rank+'/'+clade+'/regression_'+clade+'_'+feature.replace(':','').replace(' ','_')+'.png')		
		plt.cla()
		plt.clf()
		plt.close()	
		
	
		#record the values for each clade-feature-OGTs to be available for further analysis/replot
		g = open('./files/individual_features/'+rank+'/'+clade+'/raw_data_'+clade+'_'+feature.replace(':','').replace(' ','_')+'.txt','w')
		g.write('OGT\tfeature_value\n')
		for x in range(0,len(x_values),1):
			g.write(str(x_values[x])+'\t'+str(y_values[x])+'\n')
		g.close()
	return (feature,feature_r)


def rs(features,species_features,species_OGTs,rank,clade,valid_species,unit):
	#create a list of values needed for analysis
	feature_list = [(feature,species_features,species_OGTs,rank,clade,valid_species,unit) for feature in features]
	'''
	#run in parallel
	p = mp.Pool()
	results = p.map(r_calc,feature_list)
	p.close()
	p.join()
	'''
	#run single thread if getting glyph error
	results = list(map(r_calc,feature_list))
	#write out the feature correlation
	feature_rs ={x[0]:x[1] for x in results if not(x[1] == None)}
	sorted_rs = reversed(sorted(feature_rs, key=lambda dict_key: abs(feature_rs[dict_key])))
	gg = open('./files/'+rank+'_'+clade+'_feature_correlation.txt','w')
	for x in sorted_rs:
		gg.write(x)
		gg.write('\t')
		gg.write(str(feature_rs[x]))
		gg.write('\n')
	gg.close()
	return feature_rs
