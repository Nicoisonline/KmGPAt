import csv
import urllib.request

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
		
def bacterie_download_data(taxname = None, projectID = None):
	"""Return the correct format to download the file from ftp given the taxname or projectID"""
	if taxname :
		with open("data/summary.txt", "r") as f:
			reader = csv.reader(f, delimiter="\t")
			for row in reader:
				if taxname in row:
					taxname_underscore = taxname.replace(" ", "_")
					return taxname_underscore + "_uid" + row[4], row[0].split(".")[0]
		print("Taxname not found in summary.txt")
	if projectID:
		with open("data/summary.txt", "r") as f:
			reader = csv.reader(f, delimiter="\t")
			for row in reader:
				if projectID in row:
					taxname = row[5]
					taxname_underscore = taxname.replace(" ", "_")
					return taxname_underscore + "_uid" + row[4], row[0].split(".")[0]
		print("Project ID not found in summary.txt")
		
def download_bacteria(taxname = None, projectID = None):
	"""Download the file from ftp given the taxname or projectID"""
	bacteria_dl_data = bacterie_download_data(taxname, projectID)
	try:
		download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/" + bacteria_dl_data[0] + "/" + bacteria_dl_data[1] + ".faa", "data/protseq/" + bacteria_dl_data[1] + ".faa")
		print("Protseq file downloaded successfully")
	except:
		print("Error downloading protseq file : Bacteria " + taxname)
	try:
		download_file("ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/" + bacteria_dl_data[0] + "/" + bacteria_dl_data[1] + ".fna", "data/genomes/" + bacteria_dl_data[1] + ".fna")
		print("Genome file downloaded successfully")
	except:
		print("Error downloading genomes file : Bacteria " + taxname)