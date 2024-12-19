import matplotlib.pyplot as plt

def get_genome(genome_id):
	# We import the file and read it
	genome = ""
	with open("data/genomes/" + genome_id + ".fna") as f:
		genome = f.read()
	# We remove the first line
	genome = genome.split("\n")[1:]
	# We join the lines
	genome = "".join(genome)
	return genome


def create_dictionary(k):
	# We create the dictionnary with all possibles k-mers A,C,G,T as keys in alphabetical order
  # example for k = 3, dict.keys = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
  
	d = {}
	for i in range(4**k):
		kmer = ""
		for j in range(k):
			kmer += "ACGT"[i // 4**(k-j-1) % 4]
		d[kmer] = 0
	return d

def get_kmers(genome, k):
	# We create the dictionnary
	kmer_dict = create_dictionary(k)
	# For each kmers of size k in the genome
	for i in range(len(genome) - k + 1):
		kmer = genome[i:i+k]
		# We only take into account the kmers with A,C,G,T
		if "N" not in kmer:
			kmer_dict[kmer] += 1
	return kmer_dict

def show_kmers(kmers):
	# We draw with plt the histogram of the kmers
	plt.bar(range(len(kmers)), kmers.values(), align='center')
	plt.xticks(range(len(kmers)), kmers.keys())
	plt.show()

def compare_kmers_graph(genome1, genome2, k):
	# We get the kmers of size k for each genome
	kmers1 = get_kmers(genome1, k)
	kmers2 = get_kmers(genome2, k)
	# We draw the histogram of the kmers
	plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
	plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)
	plt.xticks(range(len(kmers1)), kmers1.keys())
	# Add legends for colour
	plt.legend(["NC_010163", "NC_012483"])
	plt.show()
