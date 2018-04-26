#!/usr/bin/env python

"""
This is a program that that counts gene expression from OrthoGroup output.
Input incldudes one RSEM file and the orthogroup text file.
output file is a two column list of the orthogroup name and the total number
of mapped reads for all the genes within that orthogroup.
Works best on super computer because it is slow
Created by Joanna S. Griffiths on April 2017
Copyright 2017 Joanna S. Griffiths. All rights reserved.

README:
-Rename genes.results file so that there are no periods in the name
-Remove header line in RSEM count file and then move everything up to the first line (so the first line isn't blank)
"""


import argparse
import os


def files():
    parser = argparse.ArgumentParser(
            description='input: RSEM genes.results file and orthogroup.txt'
        )
    parser.add_argument(
            '-count',
            required=True,
            help="enter RSEM genes.results file",
            type=str
        )
    parser.add_argument(
            '-ortho',
            required=True,
            help="enter orthogroup.txt file name",
            type=str
        )
    return parser.parse_args()


def make_ortho_list(args):
    with open(args.ortho, 'r') as infile:
        Ortho_List=[]
        for line in infile:
            new_list = line.replace('\n', '') #removes \n if they are the end of the line
            Ortho_List.append(new_list.split(" ")) #splits the line by spaces into a list and then adds this list to the Ortho_List
    return Ortho_List #returns a list where each item is a list of the orthogroup and all the gene names


def make_RSEM_list(args):
    with open(args.count_1, 'r') as infile:
        protein_count1 = []
        for line in infile:
            protein_list = []
            keep = line.split('\t') #split lines by tab
            protein_list.append(keep[0]) #keep the gene_id (column 1)
            protein_list.append(float(keep[4])) #keep the expected count for that gene_id and turn it into a decimal number; assuming count is in column 4
            protein_count1.append(protein_list) #add gene_id and its map count to the protein count list as a list
        return protein_count1 #returns a list where each item is a list of the protein name and its expected count as a decimal


def match_RSEM_pop1(ortholist, RSEM_list):
    mylist2=[]
    for orthogroup in ortholist: #for each list(orthogroup name and protein names) in the orthogroup list
        mylist=[]
        for group in orthogroup: #for each protein name in the orthogroup list
            for protein in RSEM_pop1: #for each list of protein and its expected count
                if protein[0] == group: #if the protein name matches the protein name in the orthogroup list, then add the expected count to a list
                    mylist.append(protein[1])
        mylist2.append(mylist) #all the proteins from one orthogroup will be in their own list, which is then added as an item to this list
    return mylist2


def sum_all_proteins_in_orthogroup(matches):
    total_count_per_orthogroup=[]
    for counts_in_othrogroup in matches: #for all the counts for each orthogroup, sum the counts, so you get a total count for each orthogroup
        total_count_per_orthogroup.append(sum(counts_in_othrogroup))
    return total_count_per_orthogroup


def make_matrix(args, match_sums, outfile_name):
    with open(args.ortho, 'r') as infile:
        Ortho_List=[]
        for line in infile:
            Ortho_List.append(line.split(None, 1)[0]) #the first word in Orthogroup added to the Ortho_List
    with open(outfile_name, 'w') as outfile:
        combined_list = [Ortho_List, match_sums]
        for x in zip(*combined_list): #prints the orthogroup name in column 1, RSEM1 counts in column 2
            outfile.write("{0}\t{1}\n".format(*x))


def main():
    args = files()
    ortholist = make_ortho_list(args) #makes the orthogroup.txt into an iterable list
    RSEM_list = make_RSEM_list(args) #makes a list of the gene_ids and its associated map counts
    matches = match_RSEM_pop1(ortholist, RSEM_list) #match the gene_id and its map count to its orthogroup
    match_sums = sum_all_proteins_in_orthogroup(matches) #sum all the counts from each gene_id for a total count for each orthogroup
    outfile_name = os.path.basename(args.count)
    outfile_name = os.path.join(outfile_name + "_orthogroup_count_matrix") #makes file name. if RSEM file is 11_gene_results_SE then outfile is 11_gene_results_SE_orthogroup_count_matrix
    make_matrix(args, match_sums, outfile_name)


if __name__ == '__main__':
    main()
