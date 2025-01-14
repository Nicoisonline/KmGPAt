# KmGPAt - Kmer Genome Proteome Analysis tool

KmGPAt aims to create a software enabling genome/proteome analysis using kmer.

## Prerequisites

- Build in Python 3.11.9
- Python Libraries :
    - matplotlib
    - numpy
    - customtkinker (requires tkinker)
    - PIL
    - pandas
    - sklearn
    - csv
    - os
    - sys
    - urllib
    - re
    - threading
    - webbrowser

# Launching the software

To launch the software, you must run the `KmGPAt_App.py` file.

## Using the software

The software offers different functionalities to analyze genomes and proteomes.
- A download tool to download data from the NCBI database.
- A choice between compare and single analysis.
    - Compare analysis allows you to compare two genomes or proteomes.
        - You can switch the representation by clicking on the "Visual Representation" button. You will have a comparison by percentage or a visual representation of the comparison.
    - Single analysis allows you to analyze a single genome or proteome.
        - You can obtain the kmer signature.
        - You can obtain an alanyse by windows. You can choose the size of the number of windows. You can also choose to save a particular window. You can also run a PCA on the windows.
            - You can change between heatmap or variance representation.
 


## Using ncbi_interactions.py as a script

You can use the script `ncbi_interactions.py` to get download the data you wish without using the software.
In the command lines you must write the ProjectId or the TaxName. You can also write "summary" to download the summary of the data.

Exemple:
"python3 ncbi_interactions.py summary" will download the summary of the data.

"python3 ncbi_interactions.py NC_000907.1" will download the data of the project with the id NC_000907.1.

"python3 ncbi_interactions.py Haemophilus influenzae Rd KW20" will download the data of the project with the taxname Haemophilus influenzae Rd KW20.