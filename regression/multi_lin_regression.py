import multiprocessing as mp
import math
import logging
logger = logging.getLogger('regression')
import sklearn
logger.info('Sklearn version: '+sklearn.__version__)
del sklearn
from sklearn import linear_model
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['axes.linewidth'] = 3
import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.metrics import r2_score

params = {
   'axes.labelsize': 25,
   'font.size': 25,
   'legend.fontsize': 20,
   'xtick.labelsize': 25,
   'ytick.labelsize': 25,
   'text.usetex': False,
   'figure.figsize': [12, 12]
   }
plt.rcParams.update(params)

#make directories for models and predictions
if not os.path.exists('./files/regression_models'):
   os.makedirs('./files/regression_models')
if not os.path.exists('./files/predictions'):
   os.makedirs('./files/predictions')

def train(inputs):
	(features,species_features,species_OGT,train_species)=inputs
	training_feature_values = []
	training_OGTs =[]
	for species in train_species:
		#test if all features in this species, if not skip this species
		if set(features).issubset(set(species_features[species].keys())):
			species_working_feature_values = []
			for feature in features:
				#create a list of feature values
				species_working_feature_values.append(species_features[species][feature])
			training_feature_values.append(species_working_feature_values)
			training_OGTs.append(species_OGT[species])				

	training_feature_values = np.array(training_feature_values,dtype='float')	
	training_OGTs = np.array(training_OGTs,dtype='float')	
	regr = linear_model.LinearRegression(fit_intercept=True, normalize=False, copy_X=True, n_jobs=1)
	regr.fit(training_feature_values, training_OGTs)
	predicts = regr.predict(training_feature_values)
	RMSE = math.sqrt(np.mean((predicts - training_OGTs) ** 2))
	r2 = r2_score(training_OGTs,predicts)
	p = float(len(features))
	n = float(training_OGTs.shape[0])
	adj_r2 = 1-((1-r2)*(n-1)/(n-p-1))
	return {'features':features,'RMSE':RMSE,'regr':regr,'adj_r2':adj_r2}

def regress(title,species_to_test,analysis,features,species_features,species_OGTs,train_test_valid,rvalues,unit):
	#build list of features to use, if it is in the analysis type, has a calculated r-value, and the r-value is >= 0.3
	all_features = {feature:rvalues[feature] for feature in features if (feature.split()[0] in analysis) and (feature in rvalues.keys()) and (abs(rvalues[feature])>=0.3)}

	#only run regression if there are correlated features for a regression
	if len(all_features)>0:
		#find the most correlated feature
		seed_feature = list(all_features.items())[0]
		for feature in list(all_features.items())[1:]:
			if abs(feature[1])>abs(seed_feature[1]): #the most positively correlated feature
				seed_feature = feature
		seed_feature = seed_feature[0]
		all_features = list(all_features.keys())
		
		#find species which fit within the given clade. generate lists of train, test, and validation species.
		all_species = [species for species in species_features.keys() if species in species_to_test]
		training_species = [species for species in all_species if train_test_valid[species] =='train']
		testing_species = [species for species in all_species if train_test_valid[species] == 'test']
				
		#only run if there are enough training and validation species to fit a line
		if (len(training_species) > 2):
			#run a first round with the single most correlated feature
			best_return = train(([seed_feature],species_features,species_OGTs,training_species))
			current_features = best_return['features']
			current_adj_r2 = best_return['adj_r2']
			current_RMSE = best_return['RMSE']
			current_model = best_return['regr']
		
			#dummy value for the first round of the run
			current_adj_r2 = 0

			#while the model improves, keep adding features
			#run as long as model improves (via r value), it is not overdetermined, and there are features to add
			while (best_return['adj_r2'] > current_adj_r2) and (len(current_features)+1 < len(all_features)) and (len(current_features)+1 < len(training_species)):
				#replace current values with new bests
				current_features = best_return['features']
				current_adj_r2 = best_return['adj_r2']
				current_RMSE = best_return['RMSE']
				current_model = best_return['regr']

				#build list of possible feature combinations, starting from previous best
				test_features = []
				for feature in all_features:
					if not(feature in current_features):
						test_features.append((current_features+[feature],species_features,species_OGTs,training_species))

				#calculate regressions in parallel
				p = mp.Pool()
				results = p.map(train,test_features)
				p.close()
				p.join()
				'''
				#single thread for testing
				results = list(map(train,test_features))
				'''
				#find the best new regression
				best_return = results[0]
				for x in results[1:]:
					if x['adj_r2'] > best_return['adj_r2']:
						best_return = x

			#catch in case exiting loop for reason other than unimproved regression
			if best_return['adj_r2']> current_adj_r2:
				current_features = best_return['features']
				current_adj_r2 = best_return['adj_r2']
				current_RMSE = best_return['RMSE']
				current_model = best_return['regr']			

			coefs ={}
			g = open('./files/regression_models/'+title+'.txt','w')
			for x in range(0,len(current_features),1):
				coefs[current_features[x]] = current_model.coef_[x]
				g.write(current_features[x]+'\t'+str(current_model.coef_[x])+'\n')
			g.write('intercept\t'+str(current_model.intercept_)+'\n')
			coefs['intercept']=current_model.intercept_
			g.close()
			
			#create a list of coefficients
			regress_features = coefs.keys()
			regress_features_wout_intercept = [feature for feature in regress_features if not(feature == 'intercept')]

			#calculate the predictions for all species

			training_reported_OGTs = []
			testing_reported_OGTs =[]
			training_pred_OGTs =[]
			testing_pred_OGTs =[]
			h = open('./files/predictions/train_species_prediction_'+title+'.txt','w')
			i = open('./files/predictions/test_species_prediction_'+title+'.txt','w')
			h.write('species\tmeasured_value\tpredicted_value\n')
			i.write('species\tmeasured_value\tpredicted_value\n')
			for species in all_species:
				#test if all features in this species, if not skip this species
				if set(regress_features_wout_intercept).issubset(set(species_features[species].keys())):
					#calculate predicted OGT
					pred_OGT = []
					for feature in regress_features_wout_intercept:
						pred_OGT.append(coefs[feature]*species_features[species][feature])
					pred_OGT.append(coefs['intercept'])
					pred_OGT = sum(pred_OGT)					
					#record predicted OGT to appropriate train/test set
					if species in testing_species:
						testing_pred_OGTs.append(pred_OGT)
						testing_reported_OGTs.append(species_OGTs[species])
						i.write(species+'\t'+str(species_OGTs[species])+'\t'+str(pred_OGT)+'\n')
					else:
						training_pred_OGTs.append(pred_OGT)
						training_reported_OGTs.append(species_OGTs[species])
						h.write(species+'\t'+str(species_OGTs[species])+'\t'+str(pred_OGT)+'\n')
			h.close()
			i.close()

			#count number of species here to ensure the count reflects only those species with all the features for the model
			logger.info('The number of training species: '+str(len(training_pred_OGTs)))
			logger.info('The number of testing species: '+str(len(testing_pred_OGTs)))

			#find bounds of plot automatically, add ~10% white space to either side
			all_points = training_reported_OGTs+testing_reported_OGTs+training_pred_OGTs+testing_pred_OGTs
			whole_range = max(all_points)-min(all_points)
			min_bound = int(min(all_points)-whole_range/10)-1 
			max_bound = int(max(all_points)+whole_range/10)+1

			#calculate final statistics
			training_pred_OGTs = np.array(training_pred_OGTs,dtype='float')	
			training_reported_OGTs = np.array(training_reported_OGTs,dtype='float')	
			testing_pred_OGTs = np.array(testing_pred_OGTs,dtype='float')	
			testing_reported_OGTs = np.array(testing_reported_OGTs,dtype='float')	

			test_RMSE = math.sqrt(np.mean((testing_pred_OGTs - testing_reported_OGTs) ** 2))
			test_r2 = r2_score(testing_reported_OGTs,testing_pred_OGTs)

			#found the best solution, log and write out coeffs
			logger.info('the best model the RMSE='+str(test_RMSE)+', R2= '+str(test_r2))

			#plot predictions
			plt.plot(range(int(min_bound-10),int(max_bound+10),1),range(int(min_bound-10),int(max_bound+10),1),'--',linewidth=3,color='#566573',alpha=1)
			plt.plot(training_reported_OGTs, training_pred_OGTs,'.',markersize=18,color='#650387',alpha=1)
			plt.plot(testing_reported_OGTs, testing_pred_OGTs,'.',markersize=18,color='#258703',alpha=1)		
			plt.xlabel('Measured '+unit,fontsize=25)
			plt.ylabel('Predicted '+unit,fontsize=25)
			plt.xlim(min_bound,max_bound)
			plt.ylim(min_bound,max_bound)
			plt.tick_params(axis='x',which='both',bottom=False,top=False)
			plt.tick_params(axis='y',which='both',left=False,right=False)
			fig_title = title+'\nN='+str(len(training_reported_OGTs)+len(testing_reported_OGTs))+' R2= '+str(test_r2)[0:6]+' RMSE='+str(test_RMSE)[0:6]
			plt.title(fig_title)		
			plt.savefig('./figures/'+title+'_regression.png')
			plt.cla()
			plt.clf()
			plt.close()
