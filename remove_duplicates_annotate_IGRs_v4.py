#!/usr/bin/env python3

'''
removeduplicates_annontateIGRS_v4.py
UPDATE: v4 2/21/23 to change the format to gff3 and updated the code to use pandas instead of csv parsing
UPDATE: v3 12/13/22, updated to allow plasmids
UPDATES: v2 1/13/21
This code removes duplicate entries from a CSV genome file and annonates all the IGRs in between genes.
How to prepare for running this code:


Compile a large .txt file of all the gff3 from the NCBI FTP website. This can be done with the conda NCBI datasets download, where you can
select the gff3 files and concatentate them into a big .txt file.

This is your input file, as a .txt. If it is provided as a .csv you will need to change the delimitor to "," in the first open statement.
The output file is in the following format: [ 'Genome','Gene', 'Symbol', 'Feature Type', 'Start', 'End', 'Strand']. There will be a new header for each new species annotated.
'''

import csv
import sys
import pandas as pd


def removeduplicates(FeatureTables, outputfile):

    with open(FeatureTables, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter='\t')

        db_genomename = []
        db_featuretype = []
        db_genestrand = []
        db_genename = []
        db_genestart =[]
        db_geneend = []
        db_locus_tag = []
        db_product = []

        for entry in fhcsv:
            ##There are some random entries the top of the gff3 file but they start with #, so we are going to ignore those
            if "#" not in entry[0]:
                    if entry[0]:

                        featuretype = entry[2]  # gene, CDS, tRNA, etc
                        genomeassembly = entry[0]  # name of organism, need to match this to the genomename from the .tab output
                        genestart = entry[3]
                        geneend = entry[4]
                        genestrand = entry[6]

                        #In the gff3 file, the final line is this whole crazy thing we need to parse out.
                        extended_info = entry[8]
                        extended_info_list = []
                        extended_info_dict = {}

                        extended_info_list = extended_info.split(";")
                        for entry in extended_info_list:
                            entry.lstrip("'").rstrip("'")
                            placeholderlist = []
                            placeholderlist = entry.split("=")
                            k = placeholderlist[0]
                            v = placeholderlist[1]
                            extended_info_dict[k] = v

                        if "locus_tag" in extended_info_dict.keys():
                            locus_tag = extended_info_dict.get("locus_tag")
                        else:
                            locus_tag = "none"

                        if "gene" in extended_info_dict.keys():
                            genename = extended_info_dict.get("gene")
                        else:
                            genename = "none"

                        if "product" in extended_info_dict.keys():
                            gene_product = extended_info_dict.get("product")
                        else:
                            gene_product = "none"

                        db_genename.append(genename)
                        db_genomename.append(genomeassembly)
                        db_featuretype.append(featuretype)
                        db_genestrand.append(genestrand)
                        db_genestart.append(genestart)
                        db_geneend.append(geneend)
                        db_product.append(gene_product)
                        db_locus_tag.append(locus_tag)

    data = [db_genomename, db_genename, db_featuretype, db_genestrand, db_genestart, db_geneend, db_locus_tag, db_product]
    df = pd.DataFrame(data) #Put it all into a pandas dataframe
    df = df.transpose()
    df.columns = (['Genome', "Gene", "FeatureType", "Strand", "Start", "End", "Locus_tag", "Product"])
    pd_db = df.drop_duplicates(['Genome', 'Start', 'End'], keep='last') #Sort out the duplicates from the pandas dataframe based on the Genome name, Start, and End
    pd_db = pd_db[pd_db.FeatureType != "region"] #We want to remove all the "regions" from the table, because those are just the entire genome (like 0-2,000,000) and not an actual feature.
    pd_db['Start'] = pd_db['Start'].astype(int)
    pd_db = pd_db.sort_values(by=["Genome", "Start"])

    #pd_db.to_csv("pd_db.csv", index=False)

    genomes = pd_db['Genome'].tolist()
    genome_set = set(genomes)
    genome_set_length = len(genome_set)
    counter = 0

    IGR_df = pd.DataFrame()
    Final_df = pd.DataFrame()

    for item in genome_set:
        counter += 1
        print("Parsing genome %s. %d/%d gff3 entires complete."%(item, counter, genome_set_length))

        sorted = pd_db[pd_db["Genome"]==item]
        sorted_length = len(sorted)

        for value in range (1, sorted_length):
            current_row = sorted.iloc[value]
            #print(current_row)

            current_row_df = pd.DataFrame({'Genome':[sorted.iloc[value, 0]], "Gene":[sorted.iloc[value, 1]], "FeatureType":[sorted.iloc[value, 2]], "Strand":[sorted.iloc[value, 3]], "Start":[sorted.iloc[value, 4]],"End":[sorted.iloc[value, 5]], "Locus_tag":[sorted.iloc[value, 6]], "Product":[sorted.iloc[value, 7]]})
            #print(current_row_df)

            current_start = sorted.iloc[value, 4]
            current_start = int(current_start)
            current_end = sorted.iloc[value, 5]
            current_end = int(current_end)
            current_genome = sorted.iloc[value, 0]

            #We are just going to go through and look at every entry BEFORE the current gene, because the gene that is after that current gene will be "before" the next gene.
            before = int(value-1)
            #before_row_df = pd.DataFrame({'Genome':[sorted.iloc[before, 0]], "Gene":[sorted.iloc[before, 1]], "FeatureType":[sorted.iloc[before, 2]], "Strand":[sorted.iloc[before, 3]], "Start":[sorted.iloc[before, 4]],"End":[sorted.iloc[before, 5]], "Locus_tag":[sorted.iloc[before, 6]], "Product":[sorted.iloc[before, 7]]})

            before_end = sorted.iloc[before, 5]
            before_end = int(before_end)

            if int(current_start) - int(before_end) > 0: #to make sure that there is in fact an empty IGR in between the gene_start
                IGR_start = int(before_end + 1)
                IGR_end = int(current_start - 1)
                if IGR_end - IGR_start > 0: #This fixes is a weird bug where the IGR end is sometimes occuring before the IGR start. You also don't want an IGR with a 0 length.
                    new_row = pd.DataFrame({'Genome':[current_genome], "Gene":["IGR"], "FeatureType":["IGR"], "Strand":["IGR"], "Start":[IGR_start],"End":[IGR_end], "Locus_tag":["none"], "Product":["none"]})
                    IGR_df = pd.concat([IGR_df, new_row], join="outer")

    #Now that we have a pandas df of all our IGRs, we want to add that back to the original dataframe and then sort it based on Genome and Start location for easy viewing.
    Final_df = pd.concat([pd_db, IGR_df], join="outer")
    Final_df['Start'] = Final_df['Start'].astype(int)
    Final_df_sorted = Final_df.sort_values(by=["Genome", "Start"])

    #I want to make a little check just to be sure that there are no double IGRs, because the pattern should always be "CDS/IGR/CDS" or "CDS/CDS" but never "IGR/IGR"
    Final_df_length = len(Final_df_sorted)
    for value in range (0, Final_df_length):
        Genome_current = Final_df_sorted.iloc[value, 0]
        FeatureType_current = Final_df_sorted.iloc[value, 2]
        before = int(value-1)
        FeatureType_before = Final_df_sorted.iloc[before, 2]

        if FeatureType_current =="IGR" and FeatureType_before =="IGR":
            print("WARNING! TWO IGRs DETECTED NEXT TO EACH OTHER! THIS SHOULD NOT HAPPEN. THE OFFENDING GENOME IS %s"%(Genome_current))

    #Export it all to csv!
    Final_df_sorted.to_csv(outputfile, index=False)


if __name__ == '__main__':
    if len(sys.argv) == 3:
         removeduplicates(sys.argv[1], sys.argv[2])
    else:
         print("Usage: 1) Input .txt file of concatentated gff3's; 2) Output file name CSV ")
         sys.exit(0)
