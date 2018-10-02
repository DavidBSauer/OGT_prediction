import logging
logger = logging.getLogger('feature_calculation')

def analysis(data,t_size):
	#calculate all features
	results ={}
	results['Density']=density(data,t_size)
	results['Length']=length(data)
	results['Coding_Noncoding']=coding_noncoding(data,t_size)
	results['Coding']=coding(data,t_size)
	results['Codon']=codon(data)
	results['Start Codon']=start_codon(data)
	results['Stop Codon']=stop_codon(data)
	results['GC']=GC(data)
	results['AG']=AG(data)
	results['Dinucleotide Fraction']=dinucleotide_freq(data)
	results['Nucleotide Fraction']=nucleotide_freq(data)
	results2 = {}
	for key in results.keys():
		if isinstance(results[key],dict):
			for subkey in results[key].keys():
				results2[key+': '+subkey]=results[key][subkey]
		else: 
			results2[key] = results[key]
	return results2


def number(data):
	#calculate the number of ORFs
	count = 0.0
	for x in data:
		count = count+1
	return count

def density(data,t_size):
	#calculate the density of ORFs in the genome
	return number(data)/t_size
	
def length(data):
	#calculate average ORF length
	length = 0.0
	for x in data:
		input_seq =data[x]
		length = length+len(input_seq)
	return length/number(data)
	
def coding_noncoding(data,t_size):
	#calculate ratio of coding/noncoding 
	length =  0.0
	for x in data:
		input_seq = data[x]
		length = length+len(input_seq)
	return length/(2*t_size-length)
	
def coding(data,t_size):
	#calculate fraction of the genome that is coding
	length =0.0
	for x in data:
		input_seq = data[x]
		length = length+len(input_seq)
	return length/t_size

def nucleotide_freq(data):
	#calculate the fraction of each nucleotide in the ORFs
	#ignore Ns
	A=0.0
	G=0.0
	T=0.0
	C=0.0
	for x in data:
		input_seq = data[x]
		input_seq = input_seq.seq
		A = A+float(input_seq.count('A'))
		G = G+float(input_seq.count('G'))
		T = T+float(input_seq.count('T'))
		C = C+float(input_seq.count('C'))
	total = A+G+T+C
	return {'A':A/total,'C':C/total,'G':G/total,'T':T/total}

def codon(data):
	#calculate the fraction of each codon in the ORFs
	counts ={x:0.0 for x in ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']}
	total = 0.0
	keys = list(counts.keys())
	keys.sort()
	for x in data:
		input_seq = data[x].seq
		input_seq = [input_seq[i:i+3] for i in range(0, len(input_seq), 3)]
		for y in keys:
			for z in input_seq:
				if y == z:
					counts[y] = counts[y]+1
	for x in keys:
		total = total + counts[x]
	return {x:counts[x]/total for x in keys}	

def start_codon(data):
	#calculate the fraction of each start codon
	counts ={x:0.0 for x in ['ATG','GTG','TTG']}
	total = 0.0
	keys = counts.keys()
	for x in data:
		input_seq = data[x].seq[:3]
		for y in keys:
			if y == input_seq:
				counts[y] = counts[y]+1
	for x in keys:
		total = total + counts[x]
	return {x:counts[x]/total for x in keys}	
	
def stop_codon(data):
	#calculate the fraction of each stop codon
	counts ={x:0.0 for x in ['TAA','TAG','TGA']}
	total = 0.0
	keys = counts.keys()
	for x in data:
		input_seq = data[x].seq[-3:]
		for y in keys:
			if y == input_seq:
				counts[y] = counts[y]+1
	for x in keys:
		total = total + counts[x]
	return {x:counts[x]/total for x in keys}	


def GC(data):
	#calculate the ORF GC fraction
	#ignore Ns
	A=0.0
	G=0.0
	T=0.0
	C=0.0
	for x in data:
		input_seq = data[x]
		input_seq = input_seq.seq
		A = A+float(input_seq.count('A'))
		G = G+float(input_seq.count('G'))
		T = T+float(input_seq.count('T'))
		C = C+float(input_seq.count('C'))
	total = A+G+T+C
	return (G+C)/total

def AG(data):
	#calculate the ORF AG fraction
	#ignore Ns
	A=0.0
	G=0.0
	T=0.0
	C=0.0
	for x in data:
		input_seq = data[x]
		input_seq = input_seq.seq
		A = A+float(input_seq.count('A'))
		G = G+float(input_seq.count('G'))
		T = T+float(input_seq.count('T'))
		C = C+float(input_seq.count('C'))
	total = A+G+T+C
	return (G+A)/total


def dinucleotide_freq(data):
	#calculate the ORF dinucleotide fraction
	#ignore dinucleotides with Ns
	counts ={x:0.0 for x in ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']}
	total = 0.0
	keys = list(counts.keys())
	keys.sort()
	for x in data:
		input_seq = data[x]
		input_seq = input_seq.seq
		for y in keys:
			counts[y] = counts[y]+float(input_seq.count(y))
	for x in keys:
		total = total + counts[x]
	return {x:counts[x]/total for x in keys}
