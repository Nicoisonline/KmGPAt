import matplotlib.pyplot as plt
import numpy as np

def get_genome_id(genome_id):
	# We import the file and read it
	genome = ""
	with open("data/genomes/" + genome_id + ".fna") as f:
		genome = f.read()
	# We remove the first line
	genome = genome.split("\n")[1:]
	# We join the lines
	genome = "".join(genome)
	return genome

def get_genome(file):
	# We import the file and read it
	genome = ""
	with open("data/genomes/" + file) as f:
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
		if kmer in kmer_dict:
			kmer_dict[kmer] += 1
	# For each kmer, we modify it into it's frequency
	for kmer in kmer_dict:
		kmer_dict[kmer] = (kmer_dict[kmer] / (len(genome) - k + 1 ))* 100
	return kmer_dict

def show_kmers(kmers):
	# We draw with plt the histogram of the kmers
	plt.bar(range(len(kmers)), kmers.values(), align='center')
	plt.xticks(range(len(kmers)), kmers.keys())
	plt.show()

def compare_kmers_graph(genome1, genome2, k, show = False, save=False, comparaison_mode = 0):
	# We get the kmers of size k for each genome
	kmers1 = get_kmers(genome1, k)
	kmers2 = get_kmers(genome2, k)

	if comparaison_mode == 0:

		kmer_diff_values = []
		# We get the difference between the two genomes
		for kmer in kmers1:
			kmer_diff_values.append(kmers1[kmer] - kmers2[kmer])

		# We draw the histogram of the difference
		plt.bar(range(len(kmers1)), kmer_diff_values, align='center')

		# We draw the histogram of the kmers with transparency
		#plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
		#plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)
	
	if comparaison_mode == 1:
		# We draw the histogram of the kmers with transparency
		plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
		plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)

		plt.legend(["NC_010163", "NC_012483"])

	# Title
	plt.title("Comparison of k-mers of size " + str(k))
	# We want kmers to be displayed on the x-axis, at -90°, with smaller font
	plt.xticks(range(len(kmers1)), kmers1.keys(), rotation=-90, fontsize=8)
	# Add legends for colour
	
	if save:
		plt.savefig("output.png")

	if show:
		plt.show()
	plt.close()

def kmer_pipeline(file1 : str, file2 : str, k : int, show=False, save=False, comparaison_mode = 0):
	#plt.ioff()
	# We get the genomes
	genome1 = get_genome_id(file1)
	genome2 = get_genome_id(file2)
	# We compare the kmers
	compare_kmers_graph(genome1, genome2, k, show, save, comparaison_mode)

def sliding_window(seq, x, k):
    """seq : main seq
    x : number of windows"""

    windows = []
    window_size = len(seq) // x

    for i in range(0, len(seq), window_size):
        windows.append(seq[i:i+window_size])

    #if the last window is smaller than the rest of the windows then add the rest of the sequence to the last window
    if len(windows[-1]) < window_size/2:
        windows[-2] += windows[-1]
        windows.pop(-1)
    return windows

def get_windows(seq, x):
    """seq : main seq
    x : number of windows"""

    windows = np.array_split(np.array(list(seq)), x)
    windows = ["".join(window) for window in windows]

    return windows

def kmer_for_windows(seq, x, k):
    """seq : main seq
    x : number of windows
    k : kmer size"""

    windows = get_windows(seq, x)
    kmers = []

    for window in windows:
        kmers.append(get_kmers(window, k))

    return kmers

def windows_heatmap(seq, x, k, show=False, save=False):
	sliding = kmer_for_windows(seq, x, k)

	val = [list(sliding[t].values()) for t in range(x)]

	# we calculate the mean of each kmer
	moy = [sum([val[j][i] for j in range(x)])/x for i in range(len(val[0]))]

	heat = [[val[j][i]-moy[i] for i in range(len(moy))] for j in range(x)]

	# we plot the heatmap x is for the number of windows
	# y is for the number of kmers, heat is a list of list of the kmer values per window
	heat_array = np.array(heat)

	plt.imshow(heat_array.T, interpolation='nearest', aspect='auto')
	plt.colorbar()

	# We change the xticks to show the windows
	plt.xticks(np.arange(x), np.arange(x))
	plt.xlabel("Windows")

	# We change the yticks to show the kmers
	kmers = create_dictionary(k)
	plt.yticks(np.arange(len(kmers)), kmers)
	# We reduce the front size of the yticks
	plt.yticks(fontsize=5)
	plt.ylabel("Kmers")

	plt.title("Heatmap of kmers in windows")

	if save:
		plt.savefig("output.png")

	if show:
		plt.show()
	plt.close()

def variance(seq, x, k, ctr, show=False, save=False):

	sliding = kmer_for_windows(seq, x, k)

	val = [list(sliding[t].values()) for t in range(x)]

	comp = np.zeros((x))
	for i in range(x):
		for j in range(i, x, 1):
			tmp = []
			for v in range(len(val[0])):
				if abs(val[i][v] - val[j][v]) > ctr:
					tmp.append(1)
				else :
					tmp.append(0)
			if sum(tmp) > 1:
				comp[i] += 1
				comp[j] += 1

	# we plot the variance histogram
	plt.bar(np.arange(x), comp)
	plt.xlabel("Windows")
	plt.xticks(np.arange(x), np.arange(x))
	kmers = create_dictionary(k)

	plt.ylabel("Number of different windows")
	plt.yticks(np.arange(0, x, 1))
	plt.title("Difference between windows")

	if save:
		plt.savefig("output.png")
	if show:
		plt.show()
	plt.close()

def kmer_single(genome : str, k : int, show=False, save=False, comparaison_mode = 0):
	kmer_single = get_kmers(genome, k)
	
	plt.bar(range(len(kmer_single)), kmer_single.values(), align='center')
	plt.title("Kmers of size " + str(k))
	plt.xticks(range(len(kmer_single)), kmer_single.keys(), rotation=-90, fontsize=8)

	if save:
		plt.savefig("output.png")
	if show:
		plt.show()
	plt.close()

def kmer_pipeline(file1 : str, file2 : str, k : int, show=False, save=False, comparaison_mode = 0):
	#plt.ioff()
	# We get the genomes
	genome1 = get_genome_id(file1)
	genome2 = get_genome_id(file2)
	# We compare the kmers
	compare_kmers_graph(genome1, genome2, k, show, save, comparaison_mode)

def single_pipeline(file : str, k : int, show=False, save=False, kmer_or_window = 0, window_size = 0, heatmap_or_variance = 0, variance_threshold = 0.1):
	genome = get_genome_id(file)
	if kmer_or_window == 0:
		kmer_single(genome, k, show, save)
	elif kmer_or_window == 1:
		if heatmap_or_variance == 0:
			windows_heatmap(genome, window_size, k, show, save)
		elif heatmap_or_variance == 1:
			variance(genome, window_size, k, variance_threshold, show, save)






# We do the same pipeline and function for amino acides
def get_protseq_id(protseq_id):
	# We import the file and read it
	protseq = ""
	with open("data/protseq/" + protseq_id + ".faa") as f:
		protseq = f.read()
	# We split the file by lines
	protseq = protseq.split("\n")
	protseqs = []

	# We are going to store each protein in a list, we disregard the first element
	# Examples for each protein
	# >name_prot1
	# SEQ1
	# ...
	# >name_prot2
	# SEQ2
	# ...
	# The result will be : [SEQ1, SEQ2, ...]
	protein = ""
	for line in protseq[1:]:
		if line.startswith(">"):
			if protein:
				protseqs.append(protein)
				protein = ""
		else:
			protein += line
	if protein:
		protseqs.append(protein)

	return protseqs

def get_protseq(file):
	# We import the file and read it
	protseq = ""
	with open("data/protseq/" + file) as f:
		protseq = f.read()
	# We split the file by lines
	protseq = protseq.split("\n")
	protseqs = []

	# We are going to store each protein in a list, we disregard the first element
	# Examples for each protein
	# >name_prot1
	# SEQ1
	# ...
	# >name_prot2
	# SEQ2
	# ...
	# The result will be : [SEQ1, SEQ2, ...]
	protein = ""
	for line in protseq[1:]:
		if line.startswith(">"):
			if protein:
				protseqs.append(protein)
				protein = ""
		else:
			protein += line
	if protein:
		protseqs.append(protein)

	return protseqs

def create_dictionary_amino_acids(k):
	d = {}
	for i in range(20**k):
		kmer = ""
		for j in range(k):
			kmer += "ACDEFGHIKLMNPQRSTVWY"[i // 20**(k-j-1) % 20]

		d[kmer] = 0
	return d

def get_kmers_amino_acids(protseq_list : list, k):
	kmer_dict = create_dictionary_amino_acids(k)
	for protseq in protseq_list:
		for i in range(len(protseq) - k + 1):
			kmer = protseq[i:i+k]
			if kmer in kmer_dict:
				kmer_dict[kmer] += 1
	for kmer in kmer_dict:
		kmer_dict[kmer] = (kmer_dict[kmer] / (len(protseq_list) - k + 1 ))* 100
	return kmer_dict

def compare_kmers_graph_amino_acids(genome1, genome2, k, show = False, save=False, comparaison_mode = 0):
	kmers1 = get_kmers_amino_acids(genome1, k)
	kmers2 = get_kmers_amino_acids(genome2, k)

	if comparaison_mode == 0:

		kmer_diff_values = []
		for kmer in kmers1:
			kmer_diff_values.append(kmers1[kmer] - kmers2[kmer])

		plt.bar(range(len(kmers1)), kmer_diff_values, align='center')

	if comparaison_mode == 1:
		plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
		plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)

		plt.legend(["NC_010163", "NC_012483"])

	plt.title("Comparison of k-mers of size " + str(k))
	plt.xticks(range(len(kmers1)), kmers1.keys(), rotation=-90, fontsize=8)

	if save:
		plt.savefig("output.png")

	if show:
		plt.show()
	plt.close()

def kmer_pipeline_amino_acids(file1 : str, file2 : str, k : int, show=False, save=False, comparaison_mode = 0):
	genome1 = get_protseq_id(file1)
	genome2 = get_protseq_id(file2)
	compare_kmers_graph_amino_acids(genome1, genome2, k, show, save, comparaison_mode)

def kmer_single_amino_acids(genome : str, k : int, show=False, save=False, comparaison_mode = 0):
	kmer_single = get_kmers_amino_acids(genome, k)
	
	plt.bar(range(len(kmer_single)), kmer_single.values(), align='center')
	plt.title("Kmers of size " + str(k))
	plt.xticks(range(len(kmer_single)), kmer_single.keys(), rotation=-90, fontsize=8)

	if save:
		plt.savefig("output.png")
	plt.close()

def sliding_window_amino_acids(seq, x, k):
	"""seq : main seq
	x : number of windows"""

	windows = []
	window_size = len(seq) // x

	for i in range(0, len(seq), window_size):
		windows.append(seq[i:i+window_size])

	if len(windows[-1]) < window_size/2:
		windows[-2] += windows[-1]
		windows.pop(-1)
	return windows

def get_windows_amino_acids(seq, x):
	"""seq : main seq
	x : number of windows"""

	windows = np.array_split(np.array(list(seq)), x)
	windows = ["".join(window) for window in windows]

	return windows

def kmer_for_windows_amino_acids(seq, x, k):
	"""seq : main seq
	x : number of windows
	k : kmer size"""

	windows = get_windows_amino_acids(seq, x)
	kmers = []

	for window in windows:
		kmers.append(get_kmers_amino_acids(window, k))

	return kmers

def protseq_heatmap(protseq_list, x, k, show=False, save=False):
	sliding = kmer_for_windows_amino_acids(protseq_list, x, k)

	val = [list(sliding[t].values()) for t in range(x)]

	moy = [sum([val[j][i] for j in range(x)])/x for i in range(len(val[0]))]

	heat = [[val[j][i]-moy[i] for i in range(len(moy))] for j in range(x)]

	heat_array = np.array(heat)

	plt.imshow(heat_array.T, interpolation='nearest', aspect='auto')
	plt.colorbar()

	plt.xticks(np.arange(x), np.arange(x))
	plt.xlabel("Windows")

	kmers = create_dictionary_amino_acids(k)
	plt.yticks(np.arange(len(kmers)), kmers)
	plt.yticks(fontsize=5)
	plt.ylabel("Kmers")

	plt.title("Heatmap of kmers in windows")

	if save:
		plt.savefig("output.png")

	if show:
		plt.show()
	plt.close()

def variance_amino_acids(protseq_list, x, k, ctr, show=False, save=False):
	
	sliding = kmer_for_windows_amino_acids(protseq_list, x, k)

	val = [list(sliding[t].values()) for t in range(x)]

	comp = np.zeros((x))
	for i in range(x):
		for j in range(i, x, 1):
			tmp = []
			for v in range(len(val[0])):
				if abs(val[i][v] - val[j][v]) > ctr:
					tmp.append(1)
				else :
					tmp.append(0)
			if sum(tmp) > 1:
				comp[i] += 1
				comp[j] += 1

	plt.bar(np.arange(x), comp)
	plt.xlabel("Windows")
	plt.xticks(np.arange(x), np.arange(x))
	kmers = create_dictionary_amino_acids(k)

	plt.ylabel("Number of different windows")
	plt.yticks(np.arange(0, x, 1))
	plt.title("Difference between windows")

	if save:
		plt.savefig("output.png")
	if show:
		plt.show()
	plt.close()

def single_pipeline_amino_acids(file : str, k : int, show=False, save=False, kmer_or_window = 0, window_size = 0, heatmap_or_variance = 0, variance_threshold = 0.1):
	genome = get_protseq_id(file)
	if kmer_or_window == 0:
		kmer_single_amino_acids(genome, k, show, save)
	elif kmer_or_window == 1:
		if heatmap_or_variance == 0:
			protseq_heatmap(genome, window_size, k, show, save=save)
		elif heatmap_or_variance == 1:
			variance_amino_acids(genome, window_size, k, variance_threshold, show, save)