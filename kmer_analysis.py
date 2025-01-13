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

def compare_kmers_graph(genome1, genome2, k, show = False, save=False, comparaison_mode = 0, name1 = "", name2 = ""):
	# We get the kmers of size k for each genome
	kmers1 = get_kmers(genome1, k)
	kmers2 = get_kmers(genome2, k)

	kmer_diff_values = []

	for kmer in kmers1:
			kmer_diff_values.append(kmers1[kmer] - kmers2[kmer])

	if comparaison_mode == 0:

		# We get the difference between the two genomes
		

		# We draw the histogram of the difference
		plt.bar(range(len(kmers1)), kmer_diff_values, align='center')

		# We draw the histogram of the kmers with transparency
		#plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
		#plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)
	
	if comparaison_mode == 1:
		# We draw the histogram of the kmers with transparency
		plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
		plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)

		plt.legend([name1, name2])


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

	# We return a list of couple of kmers and their difference
	result = [(kmer, kmer_diff_values[i]) for i, kmer in enumerate(kmers1.keys())]
	# We return a string of the result with "\n" as separator as format "kmer : difference"
	return "\n".join([str(kmer) + " : " + str(diff) for kmer, diff in result])


#def kmer_pipeline(file1 : str, file2 : str, k : int, show=False, save=False, comparaison_mode = 0):
#	#plt.ioff()
#	# We get the genomes
#	genome1 = get_genome_id(file1)
#	genome2 = get_genome_id(file2)	 
#	# We compare the kmers
#	compare_kmers_graph(genome1, genome2, k, show, save, comparaison_mode)

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

	# We return the total length of the genome and the length of each windows as a string
	return "Total length of the genome : " + str(len(seq)) + "\n" + "Length of each window : " + str(len(seq)//x)

def variance(seq, x, k, ctr, show=False, save=False):

	sliding = kmer_for_windows(seq, x, k)

	val = [list(sliding[t].values()) for t in range(x)]
	moy = [sum([val[j][i] for j in range(x)])/x for i in range(len(val[0]))]

	heat = [[val[j][i]-moy[i]for i in range(len(moy))] for j in range(x)]

	plt.figure()
	for k in range(len(heat)):
		plt.plot(np.arange(len(val[k])), heat[k], label="Window " + str(k+1))

	plt.xlabel("Kmers")
	# We change the xticks to show the kmers
	plt.xticks(np.arange(len(val[0])), sliding[0].keys(), rotation=-90, fontsize=8)
	plt.ylabel("Difference with the mean")
	plt.title("Difference between windows")

	if save:
		plt.savefig("output.png")
	if show:
		plt.show()
	plt.close()

	return "Total length of the genome : " + str(len(seq)) + "\n" + "Length of each window : " + str(len(seq)//x)

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

	# We return the kmers and their frequency as a string with "\n" as separator as format "kmer : frequency"
	return "\n".join([str(kmer) + " : " + str(freq) for kmer, freq in kmer_single.items()])

def kmer_pipeline(file1 : str, file2 : str, k : int, show=False, save=False, comparaison_mode = 0):
	#plt.ioff()
	# We get the genomes
	genome1 = get_genome_id(file1)
	genome2 = get_genome_id(file2)
	# We compare the kmers
	return compare_kmers_graph(genome1, genome2, k, show, save, comparaison_mode, file1, file2)

def single_pipeline(file : str, k : int, show=False, save=False, kmer_or_window = 0, window_size = 0, heatmap_or_variance = 0, variance_threshold = 0.1):
	genome = get_genome_id(file)
	if kmer_or_window == 0:
		return kmer_single(genome, k, show, save)
	elif kmer_or_window == 1:
		if heatmap_or_variance == 0:
			return windows_heatmap(genome, window_size, k, show, save)
		elif heatmap_or_variance == 1:
			return variance(genome, window_size, k, variance_threshold, show, save)






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
	# For each kmer, we modify it into it's frequency
	total_kmer = 0
	for kmer in kmer_dict:
		total_kmer += kmer_dict[kmer]
	for kmer in kmer_dict:
		kmer_dict[kmer] = (kmer_dict[kmer] / total_kmer ) * 100
	return kmer_dict

def get_kmers_amino_acids_windows(protseq_str : str, k):
	kmer_dict = create_dictionary_amino_acids(k)
	for i in range(len(protseq_str) - k + 1):
		kmer = protseq_str[i:i+k]
		if kmer in kmer_dict:
			kmer_dict[kmer] += 1
	# For each kmer, we modify it into it's frequency
	total_kmer = 0
	for kmer in kmer_dict:
		total_kmer += kmer_dict[kmer]
	for kmer in kmer_dict:
		kmer_dict[kmer] = (kmer_dict[kmer] / total_kmer ) * 100
	return kmer_dict

def compare_kmers_graph_amino_acids(genome1, genome2, k, show = False, save=False, comparaison_mode = 0, name1 = "", name2 = ""):
	kmers1 = get_kmers_amino_acids(genome1, k)
	kmers2 = get_kmers_amino_acids(genome2, k)

	kmer_diff_values = []

	for kmer in kmers1:
			kmer_diff_values.append(kmers1[kmer] - kmers2[kmer])

	if comparaison_mode == 0:

		plt.bar(range(len(kmers1)), kmer_diff_values, align='center')

	if comparaison_mode == 1:
		plt.bar(range(len(kmers1)), kmers1.values(), align='center', color='b', alpha=0.5)
		plt.bar(range(len(kmers2)), kmers2.values(), align='center', color='r', alpha=0.5)

		plt.legend([name1, name2])

	plt.title("Comparison of k-mers of size " + str(k))
	plt.xticks(range(len(kmers1)), kmers1.keys(), rotation=-90, fontsize=8)

	if save:
		plt.savefig("output.png")

	if show:
		plt.show()
	plt.close()

	# We return a string of the result with "\n" as separator as format "kmer : difference"
	result = [(kmer, kmer_diff_values[i]) for i, kmer in enumerate(kmers1.keys())]
	return "\n".join([str(kmer) + " : " + str(diff) for kmer, diff in result])


def kmer_pipeline_amino_acids(file1 : str, file2 : str, k : int, show=False, save=False, comparaison_mode = 0):
	genome1 = get_protseq_id(file1)
	genome2 = get_protseq_id(file2)
	return compare_kmers_graph_amino_acids(genome1, genome2, k, show, save, comparaison_mode, file1, file2)

def kmer_single_amino_acids(genome : str, k : int, show=False, save=False, comparaison_mode = 0):
	kmer_single = get_kmers_amino_acids(genome, k)
	
	plt.bar(range(len(kmer_single)), kmer_single.values(), align='center')
	plt.title("Kmers of size " + str(k))
	plt.xticks(range(len(kmer_single)), kmer_single.keys(), rotation=-90, fontsize=8)

	if save:
		plt.savefig("output.png")
	plt.close()

	return "\n".join([str(kmer) + " : " + str(freq) for kmer, freq in kmer_single.items()])

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
		kmers.append(get_kmers_amino_acids_windows(window, k))

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

	# We return the total length of the genome and the length of each windows as a string
	seq = "".join(protseq_list)
	return "Total length of the genome : " + str(len(seq)) + "\n" + "Length of each window : " + str(len(seq)//x)

def variance_amino_acids(protseq_list, x, k, show=False, save=False):
	
	sliding = kmer_for_windows_amino_acids(protseq_list, x, k)

	val = [list(sliding[t].values()) for t in range(x)]
	moy = [sum([val[j][i] for j in range(x)])/x for i in range(len(val[0]))]

	heat = [[val[j][i]-moy[i]for i in range(len(moy))] for j in range(x)]

	plt.figure()
	for k in range(len(heat)):
		plt.plot(np.arange(len(val[k])), heat[k], label="Window " + str(k+1))

	plt.xlabel("Kmers")
	# We change the xticks to show the kmers
	plt.xticks(np.arange(len(val[0])), sliding[0].keys(), rotation=-90, fontsize=8)
	plt.ylabel("Difference with the mean")
	plt.title("Difference between windows")

	if save:
		plt.savefig("output.png")
	if show:
		plt.show()
	plt.close()

	seq = "".join(protseq_list)
	return "Total length of the genome : " + str(len(seq)) + "\n" + "Length of each window : " + str(len(seq)//x)

def single_pipeline_amino_acids(file : str, k : int, show=False, save=False, kmer_or_window = 0, window_size = 0, heatmap_or_variance = 0, variance_threshold = 0.1):
	genome = get_protseq_id(file)
	if kmer_or_window == 0:
		return kmer_single_amino_acids(genome, k, show, save)
	elif kmer_or_window == 1:
		if heatmap_or_variance == 0:
			return protseq_heatmap(genome, window_size, k, show, save=save)
		elif heatmap_or_variance == 1:
			return variance_amino_acids(genome, window_size, k, show, save)

def save_windows(file_input : str, number_of_windows : int, windows_number : int, genome_or_protseq : int):
	"""From a file, we save the windows of the proteins in a file as : """

	if genome_or_protseq == 0:
		# Genome
		genome = get_genome_id(file_input)
		windows = sliding_window(genome, number_of_windows, windows_number)
		with open("data/genomes/" + file_input + "_windows" + str(number_of_windows) + "_window" + str(windows_number) + ".fna", "w") as f:
			f.write(">Window" + str(windows_number) + "\n" + windows[windows_number] + "\n")
	else:
		# Protseq
		protseq = get_protseq_id(file_input)
		windows = sliding_window_amino_acids(protseq, number_of_windows, windows_number)
		
        # For now the window of interest is a list of string, we need to convert it to a string
		window_str = ""
		for i in windows[windows_number]:
			window_str += i
		with open("data/protseq/" + file_input + "_windows" + str(number_of_windows) + "_window" + str(windows_number) + ".faa", "w") as f:
			f.write(">Window" + str(windows_number) + "\n" + window_str + "\n")


# Cimetiere des fonctions

def old_variance_amino_acids(protseq_list, x, k, ctr, show=False, save=False):
	
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