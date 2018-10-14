def analysis(data):
	#calculate all rRNA features
	results ={}
	N_counts = counter(data)
	results['Nucleotide Fraction']=nucleotide_freq(N_counts)
	results['GC']=GC(N_counts)
	results2 = {}
	for key in results.keys():
		if isinstance(results[key],dict):
			for subkey in results[key].keys():
				results2[key+': '+subkey]=results[key][subkey]
		else: 
			results2[key] = results[key]
	return results2

def counter(data):
	#count the nucleotide in the genome
	#ignore Ns
	A=0.0
	G=0.0
	T=0.0
	C=0.0
	for x in data:
		input_seq = data[x].seq
		A = A+float(input_seq.count('A'))
		G = G+float(input_seq.count('G'))
		T = T+float(input_seq.count('T'))
		C = C+float(input_seq.count('C'))
	return {'A':A,'G':G,'T':T,'C':C}

def GC(N_counts):
	#calculate the GC fraction
	total = sum(N_counts.values())
	return (N_counts['G']+N_counts['C'])/total
	
def nucleotide_freq(N_counts):
	#calculate the nucleotide fraction
	total = sum(N_counts.values())
	return {'A':N_counts['A']/total,'C':N_counts['C']/total,'G':N_counts['G']/total,'T':N_counts['T']/total}