#!/usr/bin/env python

"""
This is a program that finds OrthoGroups that only have genes from one population
Created by Joanna S. Griffiths on May 2017
Copyright 2017 Joanna S. Griffiths. All rights reserved.

"""

import argparse


def files():
    parser = argparse.ArgumentParser(
            description='OrthoGroups input'
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
            Ortho_List.append(new_list) #splits the line by spaces into a list and then adds this list to the Ortho_List
    return Ortho_List


def find_orthos(Ortho_List):
    GOL_only = []
    BOD_only = []
    for OrthoGroup in Ortho_List:
        if "BOD" not in OrthoGroup:  #If no BOD genes are found, add it to the GOL only list
            GOL_only.append(OrthoGroup)
        if "GOL" not in OrthoGroup:  #If no GOL genes are found, add it to the BOD only list
            BOD_only.append(OrthoGroup)
    return GOL_only, BOD_only


def make_orpahned_list(GOL_only, BOD_only):
    GOL_gene_list = []
    BOD_gene_list = []
    for OrthoGroup in GOL_only:
        temp_list_GOL = OrthoGroup.split(" ")
        GOL_gene_list.append(temp_list_GOL[1:]) #prints everything but first item which is the Orthogroup name
    for OrthoGroup in BOD_only:
        temp_list_BOD = OrthoGroup.split(" ")
        BOD_gene_list.append(temp_list_BOD[1:])
    with open("orphaned_GOL_list", 'a') as outfile:
            for OrthoGroups in GOL_gene_list:
                for gene in OrthoGroups:
                    outfile.write("{0}\n".format(gene))
    with open("orphaned_BOD_list", 'a') as outfile:
            for OrthoGroups in BOD_gene_list:
                for gene in OrthoGroups:
                    outfile.write("{0}\n".format(gene))


def main():
    args = files()
    Ortho_List = make_ortho_list(args)
    GOL_only, BOD_only = find_orthos(Ortho_List)
    make_orpahned_list(GOL_only, BOD_only) #puts the orphaned orthogroups into a text file


if __name__ == '__main__':
    main()
