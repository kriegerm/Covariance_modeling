#!/usr/bin/env python3
'''
UPDATE: 2/5/21
makeMSA_from_CMresults_v1.py

This code takes in three inputs:
1) A CSV file of all the TaxID's and Species Names, with a header row. This should be a list of the SAME species that you did your CM search on.
An example of this format is:
TaxID	Species Name
1302	Streptococcus gordonii
1304	Streptococcus salivarius
1309	Streptococcus mutans
1310	Streptococcus sobrinus
1311	Streptococcus agalactiae
1313	Streptococcus pneumoniae
1314	Streptococcus pyogenes

2) The results from your CM search (remeber, using the SAME species in the above file) in .tab format. Don't change anything about this output.

3) The name of your MSA file that's going to have 1's for all the species in the TaxID/SpeciesName list you used as input 1 that are also represented in your CM results,
 and 0's for all the species in that are in the input 1 list but are not represented in the CM results.

'''

import csv
import sys
import re

def makeMSA(inputTaxIDandSpeciesNameCSV, inputCMresultstab, outputMSA):
##1) Open the csv with TaxID in column 1 and Species name in column 2 (with header) and create lists with species name and taxid

    SpeciesName = []
    TaxID = []
    SpeciesNameandTaxID = []

    with open(inputTaxIDandSpeciesNameCSV, 'r') as tabfh:
        tabfhcsv = csv.reader(tabfh, delimiter=",")
        next(tabfhcsv)

        for entry in tabfhcsv:
            SpeciesName.append(entry[1])
            TaxID.append(entry[0])
            SpeciesNameandTaxID.append(entry[0] + " " + entry[1])


    print("There are %d taxID/Species Name pairs." %len(SpeciesName))

##2) Open the .tab CM search results, read through it and store the entries in lists

    CMResults_targetname = []
    CMResults_seqfrom = []
    CMResults_seqto = []
    CMResults_strand = []
    CMResults_gc = []
    CMResults_bias = []
    CMResults_score = []
    CMResults_eval = []
    CMResults_fullname = []
    CMResults_speciesname = []


    with open(inputCMresultstab, 'r') as CMresults:
        CMresultsfh = csv.reader(CMresults, delimiter="\t")
        fhcsv = next(CMresultsfh)

        for entry in CMresultsfh:
            entry = str(entry)
            entry = entry.rstrip("']")
            entry = entry.lstrip("['")
            entry = re.sub(' +', ',', entry)
            entry = entry.split(",")
            entry = list(entry)

            if "#" not in entry[0]:

                CMResults_targetname.append(entry[0])
                CMResults_seqfrom.append(entry[7])
                CMResults_seqto.append(entry[8])
                CMResults_strand.append(entry[9])
                CMResults_gc.append(entry[12])
                CMResults_bias.append(entry[13])
                CMResults_score.append(entry[14])
                CMResults_eval.append(entry[15])
                CMResults_speciesname.append(entry[17] + " " + entry[18])

                #Need to do a bit more work to parse out the full name and edit it without extra spaces or commas
                fullname = str((entry[17:-1]))
                fullname = fullname.replace("'", "")
                fullname = re.sub(',', ' ', fullname)
                fullname = re.sub(' +', ' ', fullname)
                CMResults_fullname.append(fullname)

    print("There are %d results in the .tab file." %len(CMResults_fullname))


##2) Figure out what names are representated in the CM results that are in the dictionary and which are not.

    SpeciesName_inCMResultsfullname = [] #This will be our list of SpeciesNames that are found in the entries

    ##IF YOU WANT THE OUTPUT TO BE JUST SPECIES NAME OR JUST TAXID AND NOT BOTH, FIDDLE WITH THE INPUT MARKED BELOW WITH #*
    for value in range(0, len(CMResults_fullname)):
        for number in range(0, len(SpeciesName)):
            if SpeciesName[number] in CMResults_fullname[value]:
                SpeciesName_inCMResultsfullname.append(SpeciesNameandTaxID[number]) #*CHANGE THIS TO SpeciesName or TaxID instead of SpeciesNameandTaxID

    print("%d of the results were found in the taxID/Species Name list." %(len(SpeciesName_inCMResultsfullname)))

    SpeciesName_NOTinCMResultsfullname = set(SpeciesNameandTaxID) - set(SpeciesName_inCMResultsfullname) #*
    SpeciesName_NOTinCMResultsfullname = list(SpeciesName_NOTinCMResultsfullname)

    print("%d remaining taxID/Species Names with no match in the CM results." %len(SpeciesName_NOTinCMResultsfullname))
    total = len(SpeciesName_inCMResultsfullname) + len(SpeciesName_NOTinCMResultsfullname)
    print("%d total entires in taxID/Species Name found and not found (two above values added up)." %total)



##3) Make your MSA file with 1's for the entires that have the genes and 0's for the entires that do not.
    with open(outputMSA, 'a') as fh:
        for entry in SpeciesName_inCMResultsfullname:
            fh.write(">" + entry + "\n")
            fh.write("1\n\n")
        for entry in SpeciesName_NOTinCMResultsfullname:
            fh.write(">" + entry +"\n")
            fh.write("0\n\n")


if __name__ == '__main__':
    if len(sys.argv) == 4:
        makeMSA(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("INPUT: 1) CSV of Taxid's and Species Names, 2)CM Search reesults .tab, 3) MSA file output in .txt format")
        sys.exit(0)

