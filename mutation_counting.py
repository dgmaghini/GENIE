"""
Author: Dylan Maghini
Date: March 17, 2019
Take in list of genes, count combinatorial mutation presence
"""
import sys

def make_genes_dict(input_file, gene_list, gene2, gene3):
    """
    Maps genes of interest to lists of mutation counts

    This function constructs a dictionary mapping (gene1, gene2, gene3) to a
    list of integers, where indices 0-7 represent the number of samples that are
    wt/wt/wt, wt/wt/mut, wt/mut/wt, wt/mut/mut, mut/wt/wt, mut/wt/mut,
    mut/mut/wt, mut/mut/mut for gene 1 (gene of interest), gene2 (background
    gene 1), and gene 3 (background gene 3).

    Parameters:
    input_file (str): Path to input file
    gene_list (list): List of genes to check in mutation table
    gene2 (str): Name of background gene
    gene3 (str): Name of background gene

    Returns:
    dict: mapping gene name to list of integers representing mutation types

    """
    genes_dict = {} # dictionary from gene of interest to all occurrences
    gene_to_ind = {} # dictionary from gene of interst to index in mutation list

    # open mutation table
    with open(input_file, "r") as f:
        mutations = f.readline().strip().split("\t")[1:]
        gene_to_ind[gene2] = mutations.index(gene2)
        gene_to_ind[gene3] = mutations.index(gene3)

        # initialize dictionary for each gene of interest
        for gene in gene_list:
            if gene in mutations and gene2 in mutations and gene3 in mutations:
                genes_dict[(gene, gene2, gene3)] = [0,0,0,0,0,0,0,0]

            # find index of each gene of interest in mutations list
            if gene in mutations:
                gene_to_ind[gene] = mutations.index(gene)

        for line in f:
            values = line.strip().split("\t")[1:]
            values = list(map(int, values))

            # iterate through every gene of interest in dictionary
            for (geneA, geneB, geneC) in genes_dict:
                # check that all genes are screened for
                if values[gene_to_ind[geneA]] != -1 and \
                   values[gene_to_ind[geneB]] != -1 and \
                   values[gene_to_ind[gene3]] != -1:

                   # check all combinations of three possible mutants
                   # update dictionary when appropriate
                   if values[gene_to_ind[geneA]] == 1:
                       if values[gene_to_ind[geneB]] == 1:
                           if values[gene_to_ind[geneC]] == 1: #mut/mut/mut
                               genes_dict[(geneA,geneB,geneC)][7] += 1
                           else: # mut/mut/wt
                               genes_dict[(geneA,geneB,geneC)][6] += 1
                       else: # gene 1 mut, gene2 wt
                           if values[gene_to_ind[geneC]] == 1: # mut/wt/mut
                               genes_dict[(geneA,geneB,geneC)][5] += 1
                           else: # mut/wt/wt
                               genes_dict[(geneA,geneB,geneC)][4] += 1
                   else: # gene 1 wt
                       if values[gene_to_ind[geneB]] == 1: #wt/mut
                           if values[gene_to_ind[geneC]] == 1: # wt/mut/mut
                               genes_dict[(geneA,geneB,geneC)][3] += 1
                           else: #wt/mut/wt
                               genes_dict[(geneA,geneB,geneC)][2] += 1
                       else: #wt/wt
                           if values[gene_to_ind[geneC]] == 1: # wt/wt/mut
                               genes_dict[(geneA,geneB,geneC)][1] += 1
                           else: # wt/wt/wt
                               genes_dict[(geneA,geneB,geneC)][0] += 1
    return genes_dict

def main():
    # take in inputs from user
    input_file = sys.argv[1]
    input_genes = sys.argv[2] # file with one input gene per line
    gene2 = sys.argv[3] # background gene 1
    gene3 = sys.argv[4] # background gene 2
    output_name = sys.argv[5]

    # parse through input genes file and add all to list
    gene_list = []
    with open(input_genes, "r") as f:
        for line in f:
            gene_list.append(line.strip().upper())

    genes_dict = make_genes_dict(input_file, gene_list, gene2, gene3)

    # iterate through all genes of interest, output dictionary contents to file
    with open(output_name, "w") as output:
        for (geneA, geneB, geneC) in genes_dict:
            vals = "\t".join(list(map(str, genes_dict[(geneA, geneB, geneC)])))
            output.write(geneA+"\t"+geneB+"\t"+geneC+"\t"+vals+"\n")

if __name__ == "__main__":
    main()
