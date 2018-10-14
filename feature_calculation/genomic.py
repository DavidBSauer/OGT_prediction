from Bio import SeqIO

def analysis(inputs):
	#calculate all genomic features of the genome
	(genome_file,species) = inputs
	folder = '_'.join(genome_file.split('.')[:-1])
	data = SeqIO.index('./output/genomes/'+species+'/'+folder+'/'+genome_file,'fasta')
	results ={}
	N_counts = counter(data)
	results['Nucleotide Fraction']=nucleotide_freq(N_counts)
	results['Dinucleotide Fraction']=dinucleotide_freq(data)
	results['GC']=GC(N_counts)
	results['Total Size']=t_size(N_counts)
	results['J2']=j2(data)
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

def nucleotide_freq(N_counts):
	#calculate the nucleotide fraction for each nucleotide in the genome
	total = sum(N_counts.values())
	return {'A':N_counts['A']/total,'C':N_counts['C']/total,'G':N_counts['G']/total,'T':N_counts['T']/total}

def dinucleotide_freq(data):
	#calculate the dinucleotide frequencies of the genome
	#ignore dinucleotides with Ns
	counts ={x:0 for x in ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']}
	total = 0
	keys = list(counts.keys())
	keys.sort()
	for x in data:
		input_seq = data[x].seq
		for y in keys:
			counts[y] = counts[y]+float(input_seq.count(y))
	for x in keys:
		total = total + counts[x]
	return {x:counts[x]/total for x in keys}

def GC(N_counts):
	#calculate the GC fraction of the genome
	total = sum(N_counts.values())
	return (N_counts['G']+N_counts['C'])/total

def t_size(data):
	#calculate the total size of the genome
	total = sum(N_counts.values())
	return total
	
def j2(data):
	#calculate the J2 metric of the genome
	YY = 0.0
	RR = 0.0
	YR = 0.0
	RY = 0.0
	for x in data:
		input_seq = data[x].seq
		for y in ['TT','CC','TC','CT']:
			YY = YY+input_seq.count(y)
		for y in ['AA','GG','AG','GA']:
			RR = RR+input_seq.count(y)
		for y in ['TA','TG','CA','CG']:
			YR = YR+input_seq.count(y)
		for y in ['AT','AC','GT','GC']:
			RY = RY+input_seq.count(y)
	total = YY+RR+YR+RY
	YY = YY/total
	RR = RR/total
	YR = YR/total
	RY = RY/total
	return YY+RR-YR-RY
