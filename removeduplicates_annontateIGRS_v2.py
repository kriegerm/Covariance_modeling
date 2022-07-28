#!/usr/bin/env python3

'''
removeduplicates_annontateIGRS_v2.py
UPDATES: v2 1/13/21
This code removes duplicate entries from a CSV genome file and annonates all the IGRs in between genes.

How to prepare for running this code:
Compile a large .txt file of all the Feature Tables from the NCBI FTP website. This can be done by running a shell script with the FTP addresses to these feature tables, as below:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/267/845/GCF_001267845.1_ASM126784v1/GCF_001267845.1_ASM126784v1_feature_table.txt.gz
One you download and unzip all these .txt feature tables, concatenate them into one big file.

This is your input file, as a .txt. If it is provided as a .csv you will need to change the delimitor to "," in the first open statement.


The output file is in the following format: [ 'Genome','Gene', 'Symbol', 'Feature Type', 'Start', 'End', 'Strand']. There will be a new header for each new species annotated.
'''

import csv
import sys
import pandas as pd


##This is the function we will use to make the list of indices in the next code block
def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list

def removeduplicates(FeatureTables, outputfile):

    with open(FeatureTables, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter='\t')

        db_genomename = []
        db_featuretype = []
        db_genestrand = []
        db_genename = []
        db_genestart =[]
        db_geneend = []
        db_geneshortname = []

        for entry in fhcsv:
            ##There are some random entries at the bottom of the .tab file but they start with #, so we are going to ignore those
            if "#" not in entry[0]:
                if entry[4] != "plasmid":
                    featuretype = entry[0]  # gene, CDS, tRNA, etc
                    genomeassembly = entry[6]  # name of organism, need to match this to the genomename from the .tab output
                    genestart = entry[7]
                    geneend = entry[8]
                    genestrand = entry[9]
                    genename = str(entry[13])
                    geneshortname = str(entry[14])
                    db_genename.append(genename)
                    db_genomename.append(genomeassembly)
                    db_featuretype.append(featuretype)
                    db_genestrand.append(genestrand)
                    db_genestart.append(genestart)
                    db_geneend.append(geneend)
                    db_geneshortname.append(geneshortname)

    data = [db_genomename, db_genename, db_geneshortname, db_featuretype, db_genestrand, db_genestart, db_geneend]
    df = pd.DataFrame(data) #Put it all into a pandas dataframe
    df = df.transpose()
    df.columns = (['Genome', "Gene", "Symbol", "FeatureType", "Strand", "Start", "End"])

    sorted = df.drop_duplicates(['Genome', 'Start', 'End'], keep='last') #Sort out the duplicates from the pandas dataframe based on the Genome name, Start, and End

    db2_genomename = []
    db2_featuretype = []
    db2_genestrand = []
    db2_genename = []
    db2_genestart = []
    db2_geneend = []
    db2_genesymbol = []

    db2_genomename = sorted['Genome'].tolist()
    db2_featuretype = sorted['FeatureType'].tolist()
    db2_genestrand = sorted['Strand'].tolist()
    db2_genename = sorted['Gene'].tolist()
    db2_genestart = sorted['Start'].tolist()
    db2_geneend = sorted['End'].tolist()
    db2_genesymbol = sorted['Symbol'].tolist()

    setofgenomes = list(set(db2_genomename)) #unique genomes in the genomename list

    for entry in setofgenomes:
        new_genome = []
        new_gene = []
        new_start = []
        new_end = []
        new_featuretype = []
        new_genestrand = []
        new_symbol = []
        new_strand = []

        counter = 0
        indices = get_index_positions(db2_genomename, entry) #this makes a list of all the index locations where a particular genome name is occuring in the big list of genome names
        for number in indices:
            counter = counter +1

            if counter == 1 :
                if int(db2_genestart[number]) != 1: #on the first iteration through the list, if there isn't already an entry appended at the first nucleotide, add an IGR there
                    if int(db2_genestart[number]) < int(db2_genestart[number + 1]): #if the start of the first entry isn't at the end of the genome

                    ###GOTTA FIX THIS FOR GENOMES WHERE THE FIRST ENTRY IS AT THE END OF THE GENOME!

                        new_genome.append(db2_genomename[number])
                        new_gene.append("IGR")
                        new_start.append("1")
                        newend = int(db2_genestart[number])
                        new_end.append(newend -1)
                        new_featuretype.append("IGR")
                        new_genestrand.append("")
                        new_symbol.append("IGR")
                        new_strand.append("")

            if counter == len(indices):
                new_genome.append(db2_genomename[number])
                new_gene.append(db2_genename[number])
                new_start.append(db2_genestart[number])
                new_end.append(db2_geneend[number])
                new_featuretype.append(db2_featuretype[number])
                new_genestrand.append(db2_genestrand[number])
                new_symbol.append(db2_genesymbol[number])
                new_strand.append(db2_genestrand[number])
                continue

            new_genome.append(db2_genomename[number])
            new_gene.append(db2_genename[number])
            new_start.append(db2_genestart[number])
            new_end.append(db2_geneend[number])
            new_featuretype.append(db2_featuretype[number])
            new_genestrand.append(db2_genestrand[number])
            new_symbol.append(db2_genesymbol[number])
            new_strand.append(db2_genestrand[number])

            new_genome.append(db2_genomename[number])
            new_gene.append("IGR")
            newstart = int(db2_geneend[number])
            new_start.append(newstart +1)
            newend = int(db2_genestart[number +1])
            new_end.append(newend - 1)
            new_featuretype.append("IGR")
            new_genestrand.append("")
            new_symbol.append("IGR")
            new_strand.append("")

        with open(outputfile, 'a') as fh:
            writer = csv.writer(fh)
            header = [ 'Genome','Gene', 'Symbol', 'Feature Type', 'Start', 'End', 'Strand']
            writer.writerow(header)
            for value in range(0, len(new_gene)):
                saver = []
                saver.append(new_genome[value])
                saver.append(new_gene[value])
                saver.append(new_symbol[value])
                saver.append(new_featuretype[value])
                saver.append(new_start[value])
                saver.append(new_end[value])
                saver.append(new_strand[value])
                writer.writerow(saver)

if __name__ == '__main__':
    if len(sys.argv) == 3:
         removeduplicates(sys.argv[1], sys.argv[2])
    else:
         print("Usage: 1) Input .txt file of Feature Tables; 2) Output file name CSV ")
         sys.exit(0)

