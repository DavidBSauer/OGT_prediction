def analysis(data,t_size):
	#calculate all features
	results ={}
	(N_counts,N_genes) = counter(data)
	results['Density']=density(N_genes,t_size)
	results['Length']=length(N_counts,N_genes)
	results['Coding_Noncoding']=coding_noncoding(N_counts,t_size)
	results['Coding']=coding(N_counts,t_size)
	results['Nucleotide Fraction']=nucleotide_freq(N_counts)
	results['Codon']=codon(data)
	results['Start Codon']=start_codon(data)
	results['Stop Codon']=stop_codon(data)
	results['GC']=GC(N_counts)
	results['AG']=AG(N_counts)
	results['Dinucleotide Fraction']=dinucleotide_freq(data)
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
	count = 0.0
	for x in data:
		count = count+1
		input_seq = data[x].seq
		A = A+float(input_seq.count('A'))
		G = G+float(input_seq.count('G'))
		T = T+float(input_seq.count('T'))
		C = C+float(input_seq.count('C'))
	return ({'A':A,'G':G,'T':T,'C':C},count)

def density(number,t_size):
	#calculate the density of ORFs in the genome
	return number/t_size
	
def length(N_counts,N_genes):
	#calculate average ORF length
	return sum(N_counts.values())/N_genes
	
def coding_noncoding(N_counts,t_size):
	#calculate ratio of coding/noncoding 
	length = sum(N_counts.values())
	return length/(2*t_size-length)
	
def coding(N_counts,t_size):
	#calculate fraction of the genome that is coding
	length = sum(N_counts.values())
	return length/t_size

def nucleotide_freq(N_counts):
	#calculate the nucleotide fraction
	total = sum(N_counts.values())
	return {'A':N_counts['A']/total,'C':N_counts['C']/total,'G':N_counts['G']/total,'T':N_counts['T']/total}
	
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

def GC(N_counts):
	#calculate the GC fraction
	total = sum(N_counts.values())
	return (N_counts['G']+N_counts['C'])/total

def AG(N_counts):
	#calculate the GC fraction
	total = sum(N_counts.values())
	return (N_counts['G']+N_counts['A'])/total

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
