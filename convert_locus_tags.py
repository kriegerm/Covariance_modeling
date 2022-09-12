#!/usr/bin/env python3

###################
#convert_locus_tags.py Last update: 9/11/22
#This code takes in a genbank file, a list of locus tags in a csv (one locus tag per line).
#For each of the locus tags, it searches through the gb file and outputs a new csv with the original locus tag
#next to the new locus tag.
#I wrote this code for India to help her convert troublesome new UA159 locus tags to their original format.

#Useage:
#genbank file, csv of locus tags with one per line, output csv file name
###################

import sys
from Bio import SeqIO
import csv

#Defining a function to use for whole gb search
def genbank_find(genbank, locus_csv, outputcsv):

    #Write in the locus tag csv and store each value in a list
    old_locus_tags_list = []
    with open(locus_csv, 'r', encoding='utf-8-sig') as csvfile:
        for line in csvfile:
            line.replace("\n", "")
            line.rstrip().lstrip()
            old_locus_tags_list.append(line.strip())
            
    #Iterate through the list of locus tags
    for item in old_locus_tags_list:
        results_list = []
        count = 0

        #Open the gb file and read it in using SeqIO
        with open(genbank) as input_handle:
            for gb_record in SeqIO.parse(input_handle, "genbank"):

                #Look through each record in the gb
                for f in gb_record.features:
                    tag = str(f.qualifiers.get('locus_tag', []))
                    tag = tag.lstrip("['").rstrip("']'")

                    #And see if the locus tag of the entry matches the item
                    if item == tag:
                        #We only want a single entries per locus tag (becuase there are both gene and CDS entires)
                        if count == 0:
                            #So once we find one, we'll add 1 to the count and stop this loop so the script can move onto the next
                            count = count + 1
                            new_tag = list(f.qualifiers.get('old_locus_tag', []))
                            results_list.append(tag)

                            #The "new tag" list is two values, so we have to iterate through it quickly to extract each value so they don't have a weird format
                            for item in new_tag:
                                results_list.append(item)
                            #I think it's fun to print the new entries so you have an idea of what's happening while the code is running
                            print(results_list)

                            #Write your results to output
                            with open(outputcsv, 'a') as fh:
                                writer = csv.writer(fh)
                                writer.writerow(results_list)

#This little block just handles how many input variables you give it to make sure you have the correct number of arguments
if __name__ == '__main__':
    if len(sys.argv) == 4:
         genbank_find(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: gb, locus tags csv, outputcsv")
         sys.exit(0)
