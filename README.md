
  _  __           _____ _____       _   
 | |/ /          / ____|  __ \ /\  | |  
 | ' / _ __ ___ | |  __| |__) /  \ | |_ 
 |  < | '_ ` _ \| | |_ |  ___/ /\ \| __|
 | . \| | | | | | |__| | |  / ____ \ |_ 
 |_|\_\_| |_| |_|\_____|_| /_/    \_\__|
                                        
                                        
# Kmer Genome Proteome Analysis tool

KmGPAt aims to create software enabling genome/proteome analysis using kmer.

## Using ncbi_interactions.py as a script

You can use the script `ncbi_interactions.py` to get download the data you wish without using the software.
In the command lines you must write the ProjectId or the TaxName. You can also write "summary" to download the summary of the data.

Exemple:
"python3 ncbi_interactions.py summary" will download the summary of the data.
"python3 ncbi_interactions.py NC_000907.1" will download the data of the project with the id NC_000907.1.
"python3 ncbi_interactions.py Haemophilus influenzae Rd KW20" will download the data of the project with the taxname Haemophilus influenzae Rd KW20.