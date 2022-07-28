#!/usr/bin/env python3
'''
UPDATE: 2/11/21
makeMSA_from_CMresults_v2.py - update adds the CSV file to the output

This code takes in three inputs:
1) A CSV file of all the TaxID's and Species Names, with a header row. This should be a list of the SAME species that you did your CM search on.
An example of this format is:
Genome ID	Genome Name
1000562.3	Streptococcus phocae C-4
1042404.3	Streptococcus thermophilus CNCM I-1630
1054460.4	Streptococcus pseudopneumoniae IS7493
1069533.5	Streptococcus infantarius subsp. infantarius CJ18
1074052.3	Streptococcus sobrinus TCI-9
1076934.5	Streptococcus lutetiensis 033
1123298.3	Streptococcus caballi DSM 19004
1123300.3	Streptococcus devriesei DSM 19639
1123301.3	Streptococcus didelphis DSM 15616

2) The results from your CM search (remeber, using the SAME species in the above file) in .tab format. Don't change anything about this output.

Example of input:
#target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
210007.7             -         6S_model-1           -          cm        1      190  1927949  1928136      +    no    1 0.45   0.0  216.9     5e-53 !   -
699248.3             -         6S_model-1           -          cm        1      190  1778172  1778360      +    no    1 0.45   0.0  216.1   8.4e-53 !   -
1123300.3            -         6S_model-1           -          cm        1      190  1250119  1250306      +    no    1 0.46   0.0  211.8     1e-51 !   -

3) The name of your MSA file that's going to have 1's for all the species in the TaxID/SpeciesName list you used as input 1 that are also represented in your CM results,
 and 0's for all the species in that are in the input 1 list but are not represented in the CM results.

4) the name of the CSV file you want to write to that has the MSA results (both species ID and name with 1 or 0)

'''

import csv
import sys
import re

def makeMSA(inputTaxIDandSpeciesNameCSV, inputCMresultstab, outputMSA, outputMSACSV):
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

    CMResults_seqfrom = []
    CMResults_seqto = []
    CMResults_strand = []
    CMResults_gc = []
    CMResults_bias = []
    CMResults_score = []
    CMResults_eval = []
    CMResults_fullname = []
    CMResults_speciesname = []
    CMResults_genomeID = []


    with open(inputCMresultstab, 'r') as CMresults:
        CMresultsfh = csv.reader(CMresults, delimiter="\t")
        next(CMresultsfh)

        for entry in CMresultsfh:
            entry = str(entry)
            entry = entry.rstrip("']")
            entry = entry.lstrip("['")
            entry = re.sub(' +', ',', entry)
            entry = entry.split(",")
            entry = list(entry)

            if "#" not in entry[0]:
                CMResults_speciesname.append(entry[0])
                CMResults_seqfrom.append(entry[7])
                CMResults_seqto.append(entry[8])
                CMResults_strand.append(entry[9])
                CMResults_gc.append(entry[12])
                CMResults_bias.append(entry[13])
                CMResults_score.append(entry[14])
                CMResults_eval.append(entry[15])

##IGNORE THIS STUFF BUT KEEP IT IF I NEED IT LATER
                #CMResults_genomeID.append(entry[-1]) ##THIS IS WHAT I REALLY NEED
                #Need to do a bit more work to parse out the full name and edit it without extra spaces or commas
                #fullname = str((entry[17:-1]))
                #fullname = fullname.replace("'", "")
                #fullname = re.sub(',', ' ', fullname)
                #fullname = re.sub(' +', ' ', fullname)
                #CMResults_fullname.append(fullname)

    print("There are %d results in the .tab file - THERE MAY BE MULTIPLE ENTRIES FOR ONE SPECIES." %len(CMResults_speciesname))



##2) Figure out what names are representated in the CM results that are in the dictionary and which are not.

    SpeciesName_inCMResultsfullname = [] #This will be our list of SpeciesNames that are found in the entries

    ##IF YOU WANT THE OUTPUT TO BE JUST SPECIES NAME OR JUST TAXID AND NOT BOTH, FIDDLE WITH THE INPUT MARKED BELOW WITH #*
    for value in range(0, len(CMResults_speciesname)):
        for number in range(0, len(TaxID)):

            if TaxID[number] in CMResults_speciesname[value]:
                SpeciesName_inCMResultsfullname.append(TaxID[number]) #*CHANGE THIS TO SpeciesName or TaxID instead of SpeciesNameandTaxID

    SpeciesName_inCMResultsfullname = set(SpeciesName_inCMResultsfullname)

    print("%d of the results were found in the taxID/Species Name list." %(len(SpeciesName_inCMResultsfullname)))

    SpeciesName_NOTinCMResultsfullname = set(TaxID) - set(SpeciesName_inCMResultsfullname) #*
    SpeciesName_NOTinCMResultsfullname = list(set(SpeciesName_NOTinCMResultsfullname))


    print("%d remaining taxID/Species Names with no match in the CM results." %len(SpeciesName_NOTinCMResultsfullname))
    print(SpeciesName_NOTinCMResultsfullname)

    total = len(SpeciesName_inCMResultsfullname) + len(SpeciesName_NOTinCMResultsfullname)
    print("%d total entires in taxID/Species Name found and not found (two above values added up)." %total)



##3) Make your MSA file with 1's for the entires that have the genes and 0's for the entires that do not.
    with open(outputMSA, 'a') as fh:
        for entry in SpeciesName_inCMResultsfullname:
            fh.write(">" + entry + "\n")
            fh.write("1\n")
        for entry in SpeciesName_NOTinCMResultsfullname:
            fh.write(">" + entry +"\n")
            fh.write("0\n")


    with open(outputMSACSV, 'a') as fh:
        writer = csv.writer(fh)
        header = ["GenomeID", "Species", "Present/Absent"]
        writer.writerow(header)
        for entry in SpeciesName_inCMResultsfullname:
            saver = []
            count = 0

            for value in range(0, len(SpeciesName)):
                if TaxID[value] == entry:
                    saver.append(SpeciesName[value])
                    count = count + 1

            if count == 0 :
                saver.append("No species name found")

            saver.append(entry)
            saver.append("1")

            writer.writerow(saver)

        for entry in SpeciesName_NOTinCMResultsfullname:
            saver = []
            count = 0

            for value in range(0, len(SpeciesName)):
                if TaxID[value] == entry:
                    saver.append(SpeciesName[value])
                    count = count + 1

            if count == 0 :
                saver.append("No species name found")

            saver.append(entry)
            saver.append("0")

            writer.writerow(saver)

if __name__ == '__main__':
    if len(sys.argv) == 5:
        makeMSA(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("INPUT: 1) CSV of Genome IDs and Species Names, 2)CM Search reesults .tab, 3) MSA file output in .txt format, 4) MSA output in .csv format")
        sys.exit(0)

