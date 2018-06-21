from scipy.stats import pearsonr
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import os

#calculate the cross correlation of all the features correlated with OGT. give as |R| to make all positive in figure

if not os.path.exists('./files/cross_correlation'):
   os.makedirs('./files/cross_correlation')
if not os.path.exists('./figures/cross_correlation'):
   os.makedirs('./figures/cross_correlation')


def cc((feature1,feature2,species_features)): #calculate the absolute value of the pearson correlation coefficient between two features
	x_values =[]
	y_values =[]
	for species in species_features.keys():
		if (feature1 in species_features[species].keys()) and (feature2 in species_features[species].keys()):
			x_values.append(species_features[species][feature1])
			y_values.append(species_features[species][feature2])
	r_value = abs(pearsonr(x_values,y_values)[0])
	return (feature1,feature2,r_value)	

def calc(rvalues,all_features,species_features,feature_order,clade):
	corr_features = [feature for feature in rvalues.keys() if abs(rvalues[feature])>=0.3] #find correlated features (abs(pearsonr) >= 0.3)
	logging.info('the number of correlated features in the '+clade+' dataset: '+str(len(corr_features)))

	#calculate abs(r) 
	#note: this calculates A-B and B-A, though this could be sped up by simply calculating A-B and then mapping this value to B-A, it's already fast enough
	pairwise = {}
	pairs = []#generate a list of pairs to test
	for feature1 in corr_features:
		pairwise[feature1] = {}
		for feature2 in corr_features:
			pairs.append((feature1,feature2,species_features))
	
	#calculate pairwise correlation in parallel
	p = mp.Pool()
	results = p.map(cc,pairs)
	p.close()
	
	#create a 2D dictionary of absolute r-values
	for (feature1,feature2,r_value) in results:
		pairwise[feature1][feature2]=r_value
	
	#sort the feature names for neatness
	sorted_features = []
	for f_class in feature_order:
		features = [x for x in all_features if x.split()[0] == f_class and x in corr_features]
		features.sort()
		sorted_features = sorted_features+features
	sorted_dict = {sorted_features[x]:x for x in range(0,len(sorted_features),1)}
	
	#initialize a numpy array for the values for plotting
	cross_corr = np.zeros(shape=(len(corr_features),len(corr_features)))

	#record values and 
	g = open('./files/cross_correlation/'+clade+'_feature_cross_correlation_abs_r_values.csv','w')
	g.write(' \t'+'\t'.join(sorted_features)+'\n')
	for feature1 in sorted_features:
		g.write(feature1+'\t')
		for feature2 in sorted_features:
			g.write(str(pairwise[feature1][feature2])+'\t')
			cross_corr[sorted_dict[feature1]][sorted_dict[feature2]] = pairwise[feature1][feature2]
		g.write('\n')
	g.close()

	#simplify names to feature class for clarity
	sorted_features_stub = [x.split()[0] for x in sorted_features]
	#reverse to match heatmap layout
	sorted_features_stub_r = sorted_features_stub
	sorted_features_stub_r.reverse()

	#plot the CC matrix as a heatmap
	plt.figure(figsize=(20,18), dpi=300)
	ax = plt.axes()
	plt.tick_params(axis='x', which='both',left=False,right=False)
	plt.tick_params(axis='y', which='both',left=False,right=False)
	pos = np.arange(len(sorted_features))
	pos = pos[::-1]
	ax.xaxis.tick_top()
	ax.set_yticks(pos + (1.0 / 2))
	ax.set_xticks(pos + (1.0 / 2))
	ax.set_yticklabels(['']*len(sorted_features_stub))
	ax.set_xticklabels(sorted_features_stub,rotation='vertical')
	new_cross_corr = cross_corr.transpose()
	new_cross_corr = np.flipud(new_cross_corr)
	p = ax.pcolormesh(new_cross_corr,vmin=0, vmax=1)
	cbar = plt.colorbar(p)
	plt.tight_layout()
	plt.savefig('./figures/cross_correlation/'+clade+'_cross_corr_abs_r_values.png')
	plt.clf()
	plt.close()

