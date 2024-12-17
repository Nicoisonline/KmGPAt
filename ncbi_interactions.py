import csv
import urllib.request
import sys
import os

def download_file(url, filename):
	"""Download the file from the given URL and save it locally"""
	urllib.request.urlretrieve(url, filename)

def download_summary_bacteria():
	"""Download the bacteria summary from ftp"""
	try:
		download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/summary.txt", "data/summary.txt")
		print("File downloaded successfully")
	except:
		print("Error downloading file : Bacteria summary")
		
def bacterie_download_data(input : str):
	"""Return the correct format to download the file from ftp given the taxname or projectID"""
	# We check if taxname start with "NC_"
	isProjectId = 0
	if input and input[:3] == "NC_":
		isProjectId = 1
	if isProjectId == 0:
		with open("data/summary.txt", "r") as f:
			reader = csv.reader(f, delimiter="\t")
			for row in reader:
				if input in row:
					taxname_underscore = input.replace(" ", "_")
					return taxname_underscore + "_uid" + row[4], row[0].split(".")[0]
		print("Taxname not found in summary.txt")
	else:
		with open("data/summary.txt", "r") as f:
			reader = csv.reader(f, delimiter="\t")
			for row in reader:
				if input in row:
					taxname = row[5]
					taxname_underscore = taxname.replace(" ", "_")
					return taxname_underscore + "_uid" + row[4], row[0].split(".")[0]
		print("Project ID not found in summary.txt")
		
def download_bacteria(input : str):
	"""Download the file from ftp given the taxname or projectID"""
	bacteria_dl_data = bacterie_download_data(input)
	if os.path.exists("data/protseq/" + bacteria_dl_data[1] + ".faa") == True:
		print("Protseq file already downloaded")
		return 2
	if os.path.exists("data/genomes/" + bacteria_dl_data[1] + ".fna") == True:
		print("Genomes file already downloaded")
		return 2
	try:
		download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/" + bacteria_dl_data[0] + "/" + bacteria_dl_data[1] + ".faa", "data/protseq/" + bacteria_dl_data[1] + ".faa")
		print("Protseq file downloaded successfully")
	except:
		print("Error downloading protseq file : Bacteria " + input)
	try:
		download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/" + bacteria_dl_data[0] + "/" + bacteria_dl_data[1] + ".fna", "data/genomes/" + bacteria_dl_data[1] + ".fna")
		print("Genome file downloaded successfully")
	except:
		print("Error downloading genomes file : Bacteria " + input)
		
def main():
	args = sys.argv[1:]
	if os.path.exists("data") == False:
		os.mkdir("data")
		os.mkdir("data/protseq")
		os.mkdir("data/genomes")
	if os.path.exists("data/protseq") == False:
		os.mkdir("data/protseq")
	if os.path.exists("data/genomes") == False:
		os.mkdir("data/genomes")
	if os.path.exists("data/summary.txt") == False:
		download_summary_bacteria()
	if len(args) == 0:
		print("Please provide the taxname or projectID")
		return 0
	else:
		if args[0] == "summary":
			print("Downloading bacteria summary")
			download_summary_bacteria()
		else:
			print("Downloading bacteria data")
			download_bacteria(args[0])

if __name__ == "__main__":
    main()