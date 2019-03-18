"""
Author: Dylan Maghini
Date: March 17, 2019
Generate sample by gene matrix for downstream analysis
"""

import sys
import synapseclient
import re

def make_panel_dictionary(input_file):
    """
    Creates dictionary from panel to genes screened by panel

    Parameters:
    input_file (str): path to file containing panel to gene info

    Returns:
    dict: Dictionary mapping panels to set of genes

    """
    panels = {}
    with open(input_file, 'r') as f:
        header = f.readline()
        line = f.readline()

        # iterate through file to pull out relevant info
        while line:
            split = line.split('\t')
            center = split[5]
            gene = split[3]
            tf = split[7]
            tf = tf.strip()
            feat = split[6]

            # if True/False, add gene to panel set
            if (tf == 'True') and (gene != "") and (feat == 'exon'):
                if center in panels:
                    panels[center].add(gene)
                else:
                    panels[center] = set([gene])
            line = f.readline()
    return panels


# open clinical sample sheet, pull out sample ids and assays for luad
def find_tumor_ids(input_file):
    """
    Parse out sample IDs and assays for all lung adenocarcinoma tumors

    Parameters:
    input_file (str): Path to file with sample metadata

    Returns:
    dict: Dictionary mapping lung adeno sample IDs to sample panel

    """
    luad_tumors = {}
    with open(input_file, 'r') as f:
        header = f.readline()+f.readline()+f.readline()+f.readline()+\
                f.readline()
        line = f.readline()

        # parse out oncotree code, sample ID, and sequence assay
        while line:
            split = line.strip().split('\t')
            sample_id = split[1]
            oncotree_code = split[3]
            seq_assay_id = split[5]
            if oncotree_code == 'LUAD':
                luad_tumors[sample_id] = seq_assay_id
            line = f.readline()
    return luad_tumors

def pull_sample_mutations(input_file, luad_tumors):
    """
    Parse out lung adenocarcinoma samples and mutations from mutation file

    This function parses a file with extended mutation information to pull out
    all relevant mutation information. It constructs a dictionary that maps
    sample IDs to a list of lists, where each sublist includes the complete
    data for a given mutation in that sample.

    Parameters:
    input_file (str): Path to file of mutation data
    luad_tumors (dict): dictionary mapping sample IDs to assay IDs

    Returns:
    dict: Dictionary mapping LUAD sample IDs to lists of mutation data

    """
    sample_data = {}
    with open(input_file, 'r') as f:
        f.readline()
        header = f.readline()
        line = f.readline()

        # parse important mutation data from each line
        while line:
            split = line.strip().split('\t')
            sample_id = split[15]
            hugo_symbol = split[0]
            entrez_gene = split[1]
            chrom = split[4]
            start = split[5]
            stop = split[6]
            strand = split[7]
            classification = split[8]
            variant_type = split[9]
            hgvsc = split[34]
            hgvsp = split[35]
            gene = split[47]
            feature = split[48]
            symbol = split[60]
            hgnc_id = split[62]
            swissprot = split[67]
            polyphen = split[72]
            variant_class = split[94]

            # if the sample ID has already been added, add additional mutation
            if sample_id in sample_data:
                sample_data[sample_id].append([hugo_symbol, entrez_gene, chrom,\
                                               start, stop, strand, \
                                               classification, variant_type, \
                                               hgvsc, hgvsp, gene, feature, \
                                               symbol, hgnc_id, swissprot, \
                                               polyphen, variant_class])
            # if sample ID has not already been added, add sample and mutation
            elif sample_id in luad_tumors:
                sample_data[sample_id] = [[hugo_symbol, entrez_gene, chrom, \
                                           start, stop, strand, classification,\
                                           variant_type, hgvsc, hgvsp, gene, \
                                           feature, symbol, hgnc_id, \
                                           swissprot, polyphen, variant_class]]

            line = f.readline()
    return sample_data

def make_mutations_list(sample_data):
    """
    Creates list of all mutations included in dataset

    This function creates a list of all mutations that occur in the sample set,
    and filters for certain mutation types (missense, nonsense, frameshift). It
    returns a list of all genes, codon symbols, and full mutations (i.e. KRAS,
    KRASp.G12, KRASp.G12D).

    Parameters:
    sample_data (dict): Dictionary from sample to mutation lists

    Returns:
    list: List of all included mutations

    """
    mutations = []
    # iterate through all sample mutation data
    for sample_id in sample_data:
        for mutation in sample_data[sample_id]:
            symbol = mutation[0] + mutation[9]
            mut_type = mutation[6]

            # construct codon symbol (gene name and codon mutated, not including
            # amino acid change)
            ints = re.findall('\d+', symbol)
            if len(ints) >= 1:
                codonsymbol = symbol[0:symbol.find(ints[0])] + ints[0]

            # add full mutation symbol to mutations
            if (symbol not in mutations) and (mut_type == \
                    "Missense_Mutation" or mut_type == "Nonsense_Mutation" or \
                    "Frame_Shift" in mut_type):
                mutations.append(symbol)

            # add codon symbol to mutations
            if (codonsymbol not in mutations) and (mut_type == \
                    "Missense_Mutation" or mut_type == "Nonsense_Mutation" or \
                    "Frame_Shift" in mut_type):
                mutations.append(codonsymbol)

            # add gene symbol to mutations
            if mutation[0] not in mutations and (mut_type == \
                    "Missense_Mutation" or mut_type == "Nonsense_Mutation" or \
                    "Frame_Shift" in mut_type):
                mutations.append(mutation[0])
    return mutations

# make dict of all mutations that should be detected in each panel, following
# the indexing schema in the 'mutations' list
def make_panel_to_muts(panels, mutations):
    """
    Makes dictionary mapping panel to all mutations screened by panel

    Parameters:
    panels (dict): Dictionary mapping panels to all genes screened by panel
    mutations (list): complete list of all mutations in GENIE dataset

    Returns:
    dict: Dictionary mapping panel to mutations screened in panel

    """
    panel_to_muts = {}
    for panel in panels:
        panel_to_muts[panel] = [False] * len(mutations) # initalize empty list

        # look for panel genes in mutations list, add mutation if contains gene
        for gene in panels[panel]:
            for i, mutation in enumerate(mutations):
                if gene in mutation:
                    panel_to_muts[panel][i] = True
    return panel_to_muts

def make_tumor_mutations(sample_data, mutations, panel_to_muts, luad_tumors):
    """
    Constructs dictionary mapping mutation presences for each sample

    This function parses through all mutations in each sample and constructs a
    dictionary mapping sample IDs to integers representing the presence or
    absence of each mutation, where -1 represents 'not screened', 0 represents
    'not mutated', and 1 represents a mutation.

    Parameters:
    sample_data (dict): mapping sample IDs to lists of mutation data
    mutations (list): all mutations included in dataset
    panel_to_muts (dict): mapping panels to mutations screened by panel
    luad_tumors (dict): mapping sample IDs to sample assays

    Returns:
    dict: mapping sample IDs to integers representing presence of mutations

    """
    tumor_mutations = {}

    # iterate through all mutations in all samples, update mutation info
    for sample_id in sample_data:
        tumor_mutations[sample_id] = [-1]*len(mutations)
        for mutation in sample_data[sample_id]:
            symbol = mutation[0] + mutation[9]
            mut_type = mutation[6]
            ints = re.findall('\d+', symbol)
            if len(ints) >= 1:
                codonsymbol = symbol[0:symbol.find(ints[0])] + ints[0]

            # if mutation present, update gene, symbol, and codon symbol to 1
            if (mut_type == "Missense_Mutation") or (mut_type == \
                "Nonsense_Mutation") or ("Frame_Shift" in mut_type):
                tumor_mutations[sample_id][mutations.index(symbol)] = 1
                tumor_mutations[sample_id][mutations.index(mutation[0])] = 1
                tumor_mutations[sample_id][mutations.index(codonsymbol)] = 1

        # check against panel for not screened vs not mutated
        panel_mutations = panel_to_muts[luad_tumors[sample_id]]
        for i, val in enumerate(tumor_mutations[sample_id]):
            if val == -1:
                if panel_mutations[i]==True:
                    tumor_mutations[sample_id][i] = 0
    return tumor_mutations


def main():
    username = sys.argv[1]
    password = sys.argv[2]

    # use the synapse client to download relevant files
    syn = synapseclient.Synapse()
    syn.login(username, password)
    combinedbed = syn.get('syn13251251')
    mutationsextended = syn.get('syn13251247')
    clinicalpatient = syn.get('syn13251229')

    # parse all files and extract relevant information
    panels = make_panel_dictionary(combinedbed.path) # panel to gene
    luad_tumors = find_tumor_ids(clinicalpatient.path) # sample to panel
    sample_data = pull_sample_mutations(mutationsextended.path, luad_tumors)
    mutations = make_mutations_list(sample_data) # find all mutations
    panel_to_muts = make_panel_to_muts(panels, mutations) # panel to mutation
    tumor_mutations = make_tumor_mutations(sample_data, mutations, \
                                            panel_to_muts, luad_tumors)
    # parse dict containing all mutations for each tumor, output to tsv
    with open("complete_mutations_table.txt", "w") as output:
        output.write("\t".join(["Sample"] + mutations)+"\n")
        for i in tumor_mutations:
            output.write("\t".join([i]+list(map(str, tumor_mutations[i])))+"\n")

if __name__ == "__main__":
    main()
