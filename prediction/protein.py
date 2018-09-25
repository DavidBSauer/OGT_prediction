import logging
logger = logging.getLogger('prediction')

#calculate protein parameters given an itterable of nucleotide ORF sequence
AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

def analysis(inputs):
	#calculate all proteome feature for a genome
	results ={}
	(AA_counts,results['Length']) = lenght(inputs)
	results['AA']=AA_freq(AA_counts)
	results['Charged']=percent_charged(AA_counts)
	results['Polar-Uncharged']=percent_polar_uncharged(AA_counts)
	results['Thermolabile']=percent_thermolabile(AA_counts)		
	results['EK_QH']=EKQH(AA_counts)
	results['EFMR']=EFMR(AA_counts)
	results['Polar_Charged']=polar_charged(AA_counts)
	results['Polar_Hydrophobic']=polar_hydrophobic(AA_counts)
	results['ERK']=percent_ERK(AA_counts)
	results['LK_Q']=LKQ(AA_counts)
	results['IVYWREL']=percent_IVYWREL(AA_counts)
	results['IVWL']=percent_IVWL(AA_counts)
	results['KVYWREP']=percent_KVYWREP(AA_counts)
	results['GARP']=percent_GARP(AA_counts)
	results['MFILVWYERP']=percent_MFILVWYERP(AA_counts)
	results['ILVYER']=percent_ILVYER(AA_counts)
	results['MILVWYER']=percent_MILVWYER(AA_counts)
	results['MFILVYERP']=percent_MFILVYERP(AA_counts)
	results['FILVYERP']=percent_FILVYERP(AA_counts)
	results['ILVWYEHR']=percent_ILVWYEHR(AA_counts)
	results['FILVWYERP']=percent_FILVWYERP(AA_counts)
	results['ILVWYGERKP']=percent_ILVWYGERKP(AA_counts)
	results['ILVYEHR']=percent_ILVYEHR(AA_counts)
	results2 ={}
	for key in results.keys():
		if isinstance(results[key],dict):
			for subkey in results[key].keys():
				results2[key+': '+subkey]=results[key][subkey]
		else: 
			results2[key] = results[key]
	return results2

	
def AA_count(inputs):
	#calculate the counts for each AA in the proteome
	AA_counts = {x:0.0 for x in AAs}
	count = 0
	for input in inputs:
		count = count +1
		protein = inputs[input].seq.translate()
		for x in AAs:	
			AA_counts[x]=AA_counts[x]+protein.count(x)
	mean_length = sum(AA_counts.values())/count
	return (AA_counts,mean_length)
		
def AA_freq(AA_counts):
	#calculate the average AA fraction in the proteome
	total = sum(AA_counts.values())
	freqs = {x:AA_counts[x]/total for x in AA_counts.keys()}
	return freqs

def percent_charged(AA_counts):
	#calculate the fraction of the proteome that is charged
	charged = 0.0
	counts =sum(AA_counts.values())
	for x in ['D','E','K','R']:	
		charged=charged+AA_counts[x]
	return charged/counts

def percent_polar_uncharged(AA_counts):
	#calculate the fraction of the proteome that is polar-uncharged
	polar = 0.0
	counts = sum(AA_counts.values())
	for x in ['S','T','N','Q']:	
		polar=polar+AA_counts[x]
	return polar/counts

def percent_thermolabile(AA_counts):
	#calculate the fraction of the proteome that is thermolabile AAs
	HQT = 0.0
	counts= sum(AA_counts.values())
	for x in ['H','Q','T']:	
		HQT=HQT+AA_counts[x]
	return HQT/counts

def EKQH(AA_counts):
	#calculate the ratio of EK/QH in the proteome
	EK = 0.0
	QT=0.0
	for x in ['E','K']:	
		EK=EK+AA_counts[x]
	for x in ['Q','T']:
		QT=QT+AA_counts[x]
	return EK/QT

def EFMR(AA_counts):
	#calculate the fraction of EFMR in the proteome
	EFMR = 0.0
	counts= sum(AA_counts.values())
	for x in ['E','F','M','R']:	
		EFMR=EFMR+AA_counts[x]
	return EFMR/counts

def polar_charged(AA_counts):
	#calculate the ratio of polar-uncharged to polar-charged in the proteome
	polar = 0.0
	charged = 0.0
	for x in ['S','T','N','Q']:	
		polar=polar+AA_counts[x]
	for x in ['D','E','K','R','H']:	
		charged=charged+AA_counts[x]
	return polar/charged
	
def polar_hydrophobic(AA_counts):
	#calculate the ratio of polar uncharged to hydrophobic in the proteome
	polar = 0.0
	hydrophobic = 0.0
	for x in ['S','T','N','Q']:	
		polar=polar+AA_counts[x]
	for x in ['A','V','I','L','M','F','Y','W']:	
		hydrophobic=hydrophobic+AA_counts[x]
	return polar/hydrophobic

def percent_ERK(AA_counts):
	#calculate the fraction of ERK in the proteome
	ERK = 0.0
	counts= sum(AA_counts.values())
	for x in ['E','R','K']:	
		ERK=ERK+AA_counts[x]
	for x in AAs:
		counts=counts+AA_counts[x]
	return ERK/counts
	
def LKQ(AA_counts):
	#calculate the ratio of LK/Q in the proteome
	LK = 0.0
	Q = 0.0
	for x in ['L','K']:	
		LK=LK+AA_counts[x]
	Q=AA_counts['Q']
	return LK/Q

def percent_IVYWREL(AA_counts):
	#calculate the fraction of IVYWREL in the proteome
	IVYWREL = 0.0
	counts =sum(AA_counts.values())
	for x in ['I','V','Y','W','R','E','L']:
		IVYWREL=IVYWREL+AA_counts[x]
	return IVYWREL/counts

def percent_IVWL(AA_counts):
	#calculate the fraction of IVWL in the proteome
	IVWL = 0.0
	counts = sum(AA_counts.values())
	for x in ['I','V','W','L']:
		IVWL=IVWL+AA_counts[x]
	return IVWL/counts

def percent_KVYWREP(AA_counts):
	#calculate the fraction of KVYWREP in the proteome
	KVYWREP = 0.0
	counts =sum(AA_counts.values())
	for x in ['K','V','W','R','E','P']:
		KVYWREP=KVYWREP+AA_counts[x]
	return KVYWREP/counts

def percent_GARP(AA_counts):
	#calculate the fraction of GARP in the proteome
	GARP = 0.0
	counts =sum(AA_counts.values())
	for x in ['G','A','R','P']:
		GARP=GARP+AA_counts[x]
	return GARP/counts
	
def percent_MFILVWYERP(AA_counts):
	#calculate the fraction of MFILVWYERP in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['M','F','I','L','V','W','Y','E','R','P']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_ILVYER(AA_counts):
	#calculate the fraction of ILVYER in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['I','L','V','Y','E','R']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_MILVWYER(AA_counts):
	#calculate the fraction of MILVWYER in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['M','I','L','V','W','Y','E','R']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_MFILVYERP(AA_counts):
	#calculate the fraction of MFILVYERP in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['M','F','I','L','V','Y','E','R','P']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_FILVYERP(AA_counts):
	#calculate the fraction of FILVYERP in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['F','I','L','V','Y','E','R','P']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_ILVWYEHR(AA_counts):
	#calculate the fraction of ILVWYEHR in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['I','L','V','W','Y','E','H','R']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_FILVWYERP(AA_counts):
	#calculate the fraction of FILVWYERP in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['F','I','L','V','W','Y','E','R','P']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_ILVWYGERKP(AA_counts):
	#calculate the fraction of ILVWYGERKP in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['I','L','V','W','Y','G','E','R','K','P']:
		counts=counts+AA_counts[x]
	return counts/total

def percent_ILVYEHR(AA_counts):
	#calculate the fraction of ILVYEHR in the proteome
	counts = 0.0
	total =sum(AA_counts.values())
	for x in ['I','L','V','Y','E','H','R']:
		counts=counts+AA_counts[x]
	return counts/total
