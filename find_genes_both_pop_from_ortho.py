#!/usr/bin/env python

"""
This is a program that that finds OrthoGroups that only have genes from both populations
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
    Both= []
    temp_both =[]
    GOL_only = []
    BOD_only = []
    for OrthoGroup in Ortho_List:
        if "BOD" not in OrthoGroup:  #If no BOD genes are found, add it to the GOL only list
            GOL_only.append(OrthoGroup)
        else:
        #if "BOD" and "GOL" in OrthoGroup:  #If no BOD genes are found, add it to the GOL only list
            temp_both.append(OrthoGroup)
    for OrthoGroup in temp_both:
        if "GOL" not in OrthoGroup:  #If no GOL genes are found, add it to the BOD only list
            BOD_only.append(OrthoGroup)
        else:
            Both.append(OrthoGroup)
    return Both


def make_orpahned_list(Both):
    with open("both_GOLBOD_ortholist", 'a') as outfile:
            for OrthoGroups in Both:
                outfile.write("{0}\n".format(OrthoGroups))


def main():
    args = files()
    Ortho_List = make_ortho_list(args)
    Both = find_orthos(Ortho_List)
    make_orpahned_list(Both)
    #print(Ortho_List)


if __name__ == '__main__':
    main()
