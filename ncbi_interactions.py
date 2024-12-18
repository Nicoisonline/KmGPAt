import csv
import urllib.request
import sys
import os
import re

# This script is used to download the bacteria summary and the bacteria data from the ftp of NCBI

# List of return codes:
# -1: No argument provided
# 0: Download successful
# 1: Error downloading Protseq file
# 2: Error downloading Genomes file
# 12 : Error downloading both Protseq and Genomes files

error = 0

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
					motif = r"subsp\.\s+\S+ "
					input = re.sub(motif, "", input)
					taxname_underscore = input.replace(" str.", "").replace(" sp.", "").replace(" ", "_")
					return taxname_underscore + "_uid" + row[4], row[0].split(".")[0]
		print("Taxname not found in summary.txt")
	else:
		with open("data/summary.txt", "r") as f:
			reader = csv.reader(f, delimiter="\t")
			for row in reader:
				if input in row:
					taxname = row[5]
					motif = r"subsp\.\s+\S+ "
					taxname = re.sub(motif, "", taxname)
					taxname_underscore = taxname.replace(" str.", "").replace(" sp.", "").replace(" ", "_")
					return taxname_underscore + "_uid" + row[4], row[0].split(".")[0]
		print("Project ID not found in summary.txt")
		
def download_bacteria(input : str):
	"""Download the file from ftp given the taxname or projectID"""
	global error
	bacteria_dl_data = bacterie_download_data(input)
	print(bacteria_dl_data[0])
	if os.path.exists("data/protseq/" + bacteria_dl_data[1] + ".faa") == True:
		print("Protseq file already downloaded")
		error = 1
	if os.path.exists("data/genomes/" + bacteria_dl_data[1] + ".fna") == True:
		print("Genomes file already downloaded")
		if error == 1:
			error = 12
		else:
			error = 2
	if error !=1 and error != 12:
		try:
			download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/" + bacteria_dl_data[0] + "/" + bacteria_dl_data[1] + ".faa", "data/protseq/" + bacteria_dl_data[1] + ".faa")
			print("Protseq file downloaded successfully")
		except:
			print("Error downloading protseq file : Bacteria " + input)
	if error !=2 and error != 12:
		try:
			download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/" + bacteria_dl_data[0] + "/" + bacteria_dl_data[1] + ".fna", "data/genomes/" + bacteria_dl_data[1] + ".fna")
			print("Genome file downloaded successfully")
		except:
			print("Error downloading genomes file : Bacteria " + input)
	return error
		
def main():
	global error
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
		return -1
	else:
		if args[0] == "summary":
			print("Downloading bacteria summary")
			download_summary_bacteria()
		else:
			print("Downloading bacteria data")
			download_bacteria(args[0])
	return error

if __name__ == "__main__":
    main()