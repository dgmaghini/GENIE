# GENIE

### Description
These scripts parse version 4 of the GENIE database. Users input a set of oncotree codes (LUAD, LUSC, etc.) representing cancer types of interest, and `parse_genie.py` downloads the relevant files and creates a sample by mutation count matrix indicating whether a particular sample is wild-type (0), mutated (1) or not screened (-1) for that mutation. `mutation_counting.py` takes the complete mutation table output from `parse_genie.py`, a list of mutations of interest (plain text file, one gene/mutation per line) and other background mutations of interest (KRAS, EGFR, etc.). This script outputs a tab-delimited file indicating the number of samples that fall into every combination of wild type and mutated for the gene of interest and background mutations. 

### Prerequisities
* Python 3.6
* synapseclient
* pysftp
* an account at [Synapse.org](https://www.synapse.org/)

### Installing
* synapseclient is easiest installed with [conda](https://anaconda.org/bioconda/synapseclient) or can be installed with `pip install synapseclient`
* parse_genie.py and mutation_counting.py can be downloaded and run using Python3

### Running

`parse_genie.py` requires a log in to the Synapse database. 
Run `parse_genie.py` with:
`python3 parse_genie.py -u <synapse_username> -p <synapse_password> -c <oncotree code 1> <oncotree code 2> ...`

`mutation_counting.py` 
Run `mutation_counting.py` with:
`python3 mutation_counting.py path/to/complete_mutations*.txt path/to/listofgenes.txt output_name.txt -g1 background_gene1 -g2 background_gene2`

### Authors
Dylan Maghini, [(view LinkedIn)](https://www.linkedin.com/in/dylan-maghini-110139101/)

### Acknowledgments

Thank you to Monte Winslow and the members of the Winslow lab for critical feedback and ideas on how to extend the use of these scripts. 
