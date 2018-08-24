#!/usr/bin/env python

"""
This is a program that gets just the first gene names from the Orthogroups
Created by Joanna S. Griffiths on April 2017
Copyright 2017 Joanna S. Griffiths. All rights reserved.

"""

import argparse
import os



def files():
    parser = argparse.ArgumentParser(
            description='Get just gene names from Orthogroups'
        )
    parser.add_argument(
            '-ortho',
            required=True,
            help="enter orthogroup.txt file name",
            type=str
        )
    parser.add_argument(
            '-out',
            required=True,
            help="enter name for outfile",
            type=str
        )
    return parser.parse_args()


def make_ortho_list(args):
    with open(args.ortho, 'r') as infile:
        Ortho_List=[]
        for line in infile:
            new_list = line.replace('\n', '') #removes \n if they are the end of the line
            Ortho_List.append(new_list.split(" ")[1]) #splits the line by spaces into a list and then adds this list to the Ortho_List
    return Ortho_List #returns a list where each item is a list of the orthogroup and all the protein names


def make_outfile(ortholist, args):
    with open(args.out, 'w') as outfile:
        for x in ortholist: #prints the orthogroup name in column 1, RSEM1 counts in column 2
            	outfile.write("{0}\n".format(x))


def main():
    args = files()
    ortholist = make_ortho_list(args)
    make_outfile(ortholist, args)


if __name__ == '__main__':
    main()
