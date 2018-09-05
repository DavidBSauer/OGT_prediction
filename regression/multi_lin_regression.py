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
from scipy.stats import pearsonr
import os

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
	(features,species_features,species_OGT,testing,valid_species)=inputs
	training_feature_values = []
	testing_feature_values =[]
	training_OGTs =[]
	testing_OGTs =[]
	for species in valid_species:
		#test if all features in this species, if not skip this species
		if set(features).issubset(set(species_features[species].keys())):
			species_working_feature_values = []
			for feature in features:
				#create a list of feature values
				species_working_feature_values.append(species_features[species][feature])
			#assign values and OGTs to either training or testing lists as appropriate
			if species in testing:
				testing_feature_values.append(species_working_feature_values)
				testing_OGTs.append(species_OGT[species])
			else:
				training_feature_values.append(species_working_feature_values)
				training_OGTs.append(species_OGT[species])				
	testing_feature_values = np.array(testing_feature_values,dtype='float')	
	training_feature_values = np.array(training_feature_values,dtype='float')	
	testing_OGTs = np.array(testing_OGTs,dtype='float')
	training_OGTs = np.array(training_OGTs,dtype='float')	
	regr = linear_model.LinearRegression(fit_intercept=True, normalize=False, copy_X=True, n_jobs=1)
	regr.fit(training_feature_values, training_OGTs)
	predicts = regr.predict(testing_feature_values)
	RMSE = math.sqrt(np.mean((predicts - testing_OGTs) ** 2))
	r = pearsonr(predicts, testing_OGTs)[0]
	return {'features':features,'r':r,'RMSE':RMSE,'regr':regr}

def regress(title,species_to_test,analysis,features,species_features,species_OGTs,testing,rvalues,unit):

	#build list of features to use, if it is in the analysis type, has a calculated r-value, and the r-value is >= 0.3
	all_features = {feature:rvalues[feature] for feature in features if (feature.split()[0] in analysis) and (feature in rvalues.keys()) and (abs(rvalues[feature])>=0.3)}

	#only run regression if there are correlated features for a regression
	if len(all_features)>0:
		#find the most correlated feature
		seed_feature = list(all_features.items())[0]
		for feature in list(all_features.items())[1:]:
			if feature[1]>seed_feature[1]: #the most positively correlated feature
				seed_feature = feature
		seed_feature = seed_feature[0]
		all_features = list(all_features.keys())
		#find species which fit within the given clade
		valid_species = [species for species in species_features.keys() if species in species_to_test]
		testing_species = [species for species in valid_species if species in testing]
		training_species = [species for species in valid_species if not(species in testing)]
		
		#only run if there are enough training and testing species to fit a line
		if (len(training_species) > 2) and (len(testing_species) > 2):
			#run a first round with the single most correlated feature
			best_return = train(([seed_feature],species_features,species_OGTs,testing,valid_species))
			current_features = best_return['features']
			current_r = best_return['r']
			current_RMSE = best_return['RMSE']
			current_model = best_return['regr']
		
			#dummy value for the first round of the run
			current_r = 0

			#while the model improves, keep adding features
			#run as long as model improves (via r value), it is not overdetermined, and there are features to add
			while (best_return['r'] > current_r) and (len(current_features)+1 < len(all_features)) and (len(current_features)+1 < len(training_species)):
				#replace current values with new bests
				current_features = best_return['features']
				current_r = best_return['r']
				current_RMSE = best_return['RMSE']
				current_model = best_return['regr']
				
				#build list of possible feature combinations, starting from previous best
				test_features = []
				for feature in all_features:
					if not(feature in current_features):
						test_features.append((current_features+[feature],species_features,species_OGTs,testing,valid_species))

				#calculate regressions in parallel
				p = mp.Pool()
				results = p.map(train,test_features)
				p.close()

				#find the best new regression
				best_return = results[0]
				for x in results[1:]:
					if x['r'] > best_return['r']:
						best_return = x

			#catch in case exiting loop for reason other than unimproved regression
			if best_return['r']> current_r:
				current_features = best_return['features']
				current_r = best_return['r']
				current_RMSE = best_return['RMSE']
				current_model = best_return['regr']			

			#found the best solution, log and write out coeffs
			logger.info('the best model the RMSE was ='+str(current_RMSE)+', r ='+str(current_r))
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
			for species in valid_species:
				#test if all features in this species, if not skip this species
				if set(regress_features_wout_intercept).issubset(set(species_features[species].keys())):
					#calculate predicted OGT
					pred_OGT = []
					for feature in regress_features_wout_intercept:
						pred_OGT.append(coefs[feature]*species_features[species][feature])
					pred_OGT.append(coefs['intercept'])
					pred_OGT = sum(pred_OGT)					
					#record predicted OGT to appropriate train/test set
					if species in testing:
						testing_pred_OGTs.append(pred_OGT)
						testing_reported_OGTs.append(species_OGTs[species])
						i.write(species+'\t'+str(species_OGTs[species])+'\t'+str(pred_OGT)+'\n')
					else:
						training_pred_OGTs.append(pred_OGT)
						training_reported_OGTs.append(species_OGTs[species])
						h.write(species+'\t'+str(species_OGTs[species])+'\t'+str(pred_OGT)+'\n')
			h.close()
			i.close()
			
			#find bounds of plot automatically
			min_bound = int(0.95*min([min(training_reported_OGTs),min(testing_reported_OGTs),min(training_pred_OGTs),min(testing_pred_OGTs)]))-1
			max_bound = int(1.05*max([max(training_reported_OGTs),max(testing_reported_OGTs),max(training_pred_OGTs),max(testing_pred_OGTs)]))+1
			
			#plot predictions
			plt.plot(range(int(min_bound*0.9),int(max_bound*1.1),1),range(int(min_bound*0.9),int(max_bound*1.1),1),'--',linewidth=3,color='#566573',alpha=1)
			plt.plot(training_reported_OGTs, training_pred_OGTs,'.',markersize=18,color='#650387',alpha=1)
			plt.plot(testing_reported_OGTs, testing_pred_OGTs,'.',markersize=18,color='#258703',alpha=1)		
			plt.xlabel('Measured '+unit,fontsize=25)
			plt.ylabel('Predicted '+unit,fontsize=25)
			#plt.xlim(0,115)
			#plt.ylim(0,115)
			#plt.xticks([20,40,60,80,100])
			#plt.yticks([20,40,60,80,100])			
			plt.tick_params(axis='x',which='both',bottom=False,top=False)
			plt.tick_params(axis='y',which='both',left=False,right=False)
			fig_title = title+'\nN='+str(len(training_reported_OGTs)+len(testing_reported_OGTs))+' r = '+str(current_r)[0:6]+' RMSE='+str(current_RMSE)[0:6]
			plt.title(fig_title)		
			plt.savefig('./figures/'+title+'.png')
			plt.cla()
			plt.clf()
			plt.close()
