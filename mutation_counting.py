"""
Author: Dylan Maghini
Date: March 17, 2019
Take in list of genes, count combinatorial mutation presence
"""

import sys
import argparse

def make_genes_dict_triples(input_file, gene_list, gene2, gene3):
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

def make_genes_dict_doubles(input_file, gene_list, gene2):
    """
    Maps genes of interest to lists of mutation counts

    This function constructs a dictionary mapping (gene1, gene2) to a
    list of integers, where indices 0-3 represent the number of samples that are
    wt/wt, wt/mut, mut/wt, mut/mut for gene 1 (gene of interest) and gene2
    (background gene 1).

    Parameters:
    input_file (str): Path to input file
    gene_list (list): List of genes to check in mutation table
    gene2 (str): Name of background gene

    Returns:
    dict: mapping gene name to list of integers representing mutation types

    """
    genes_dict = {}
    gene_to_ind = {}
    with open(input_file, "r") as f:
        mutations = f.readline().strip().split("\t")[1:]
        gene_to_ind[gene2] = mutations.index(gene2)
        for gene in gene_list:
            if gene in mutations and gene2 in mutations:
                genes_dict[(gene, gene2)] = [0,0,0,0]
            if gene in mutations:
                gene_to_ind[gene] = mutations.index(gene)

        for line in f:
            values = line.strip().split("\t")[1:]
            values = list(map(int, values))
            for (geneA, geneB) in genes_dict:
                if values[gene_to_ind[geneA]] != -1 and \
                   values[gene_to_ind[geneB]] != -1: # if were screened
                    if values[gene_to_ind[geneA]] == 1:
                        if values[gene_to_ind[geneB]] == 1:
                            genes_dict[(geneA, geneB)][3] += 1
                        else:
                            genes_dict[(geneA, geneB)][2] += 1
                    else:
                        if values[gene_to_ind[geneB]] == 1:
                            genes_dict[(geneA, geneB)][1] += 1
                        else:
                            genes_dict[(geneA, geneB)][0] += 1
    return genes_dict

def make_genes_dict_singles(input_file, gene_list):
    """
    Maps genes of interest to lists of mutation counts

    This function constructs a dictionary mapping gene1 to a list of integers,
    where indices 0-1 represent the number of samples that are wt or mut,
    respectively for gene 1 (gene of interest).

    Parameters:
    input_file (str): Path to input file
    gene_list (list): List of genes to check in mutation table

    Returns:
    dict: mapping gene name to list of integers representing mutation types

    """
    genes_dict = {}
    gene_to_ind = {}
    with open(input_file, "r") as f:
        mutations = f.readline().strip().split("\t")[1:]
        for gene in gene_list:
            if gene in mutations:
                genes_dict[gene] = [0,0] #wt, mut
                gene_to_ind[gene] = mutations.index(gene)
        for line in f:
            values = line.strip().split("\t")[1:]
            values = list(map(int, values))
            for gene in genes_dict:
                if values[gene_to_ind[gene]] != -1:
                    if values[gene_to_ind[gene]] == 1:
                        genes_dict[gene][1] += 1
                    else:
                        genes_dict[gene][0] += 1
    return genes_dict

def main():
    # take in inputs from user
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Mutation table file name")
    parser.add_argument("input_genes", help="File containing list of genes")
    parser.add_argument("output_name", help="Output file name")
    parser.add_argument("-g1", type = str, dest = "gene2", \
                        help="Second gene of interest")
    parser.add_argument("-g2", type = str, dest = "gene3", \
                        help="Second gene of interest")
    args = parser.parse_args()

    # parse through input genes file and add all to list
    gene_list = []
    with open(args.input_genes, "r") as f:
        for line in f:
            gene_list.append(line.strip().upper())

    if args.gene2:
        if args.gene3:
            genes_dict = make_genes_dict_triples(args.input_file, gene_list, \
            args.gene2, args.gene3)
            header = "\t".join(["gene1", "gene2", "gene3", "wt/wt/wt", \
                        "wt/wt/-", "wt/-/wt", "wt/-/-", "-/wt/wt", "-/wt/-", \
                        "-/-/wt", "-/-/-"])
        else:
            genes_dict = make_genes_dict_doubles(args.input_file, gene_list, \
                                    args.gene2)
            header = "\t".join(["gene1", "gene2", "wt/wt", "wt/-", "-/wt", \
                                    "-/-"])
    else:
        genes_dict = make_genes_dict_singles(args.input_file, gene_list)
        header = "\t".join(["gene1", "wt", "-"])

    # iterate through all genes of interest, output dictionary contents to file
    with open(args.output_name, "w") as output:
        output.write(header + "\n")
        for item in genes_dict:
            if args.gene2:
                gene_names = "\t".join(item)
            else:
                gene_names = item
            vals = "\t".join(list(map(str, genes_dict[item])))
            output.write(gene_names +"\t"+vals+"\n")

if __name__ == "__main__":
    main()
