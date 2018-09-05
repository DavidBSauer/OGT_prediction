from Bio import SeqIO
import logging

def analysis(inputs):
	#calculate all genomic features of the genome
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	data = SeqIO.index('./output/genomes/'+species+'/'+folder+'/'+genome_file,'fasta')
	results ={}
	results['Nucleotide Fraction']=nucleotide_freq(data)
	results['Dinucleotide Fraction']=dinucleotide_freq(data)
	results['GC']=GC(data)
	results['Total Size']=t_size(data)
	results['J2']=j2(data)
	results2 = {}
	for key in results.keys():
		if isinstance(results[key],dict):
			for subkey in results[key].keys():
				results2[key+': '+subkey]=results[key][subkey]
		else: 
			results2[key] = results[key]
	return results2

def nucleotide_freq(data):
	#calculate the nucleotide fraction for each nucleotide in the genome
	#ignore Ns
	A=0.0
	G=0.0
	T=0.0
	C=0.0
	for x in data:
		input = data[x]
		A = A+float(input.seq.count('A'))
		G =	G+float(input.seq.count('G'))
		T = T+float(input.seq.count('T'))
		C = C+float(input.seq.count('C'))
	total = A+G+T+C
	return {'A':A/total,'C':C/total,'G':G/total,'T':T/total}

def dinucleotide_freq(data):
	#calculate the dinucleotide frequencies of the genome
	#ignore dinucleotides with Ns
	counts ={x:0 for x in ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']}
	total = 0
	keys = counts.keys()
	keys.sort()
	for x in data:
		input = data[x]
		for y in keys:
			counts[y] = counts[y]+float(input.seq.count(y))
	for x in keys:
		total = total + counts[x]
	return {x:counts[x]/total for x in keys}

def GC(data):
	#calculate the GC fraction of the genome
	#ignore Ns
	A=0.0
	G=0.0
	T=0.0
	C=0.0
	for x in data:
		input = data[x]
		A = A+float(input.seq.count('A'))
		G =	G+float(input.seq.count('G'))
		T = T+float(input.seq.count('T'))
		C = C+float(input.seq.count('C'))
	total = A+G+T+C
	return (G+C)/total

def t_size(data):
	#calculate the total size of the genome
	length = 0.0
	for x in data:
		input=data[x]
		length = length+len(input)
	return length
	
def j2(data):
	#calculate the J2 metric of the genome
	YY = 0.0
	RR = 0.0
	YR = 0.0
	RY = 0.0
	for x in data:
		input = data[x]
		for y in ['TT','CC','TC','CT']:
			YY = YY+input.seq.count(y)
		for y in ['AA','GG','AG','GA']:
			RR = RR+input.seq.count(y)
		for y in ['TA','TG','CA','CG']:
			YR = YR+input.seq.count(y)
		for y in ['AT','AC','GT','GC']:
			RY = RY+input.seq.count(y)
	total = YY+RR+YR+RY
	YY = YY/total
	RR = RR/total
	YR = YR/total
	RY = RY/total
	return YY+RR-YR-RY
