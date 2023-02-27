#!/usr/bin/env python3

'''
annonate_cmsearch_results_v4.py
UPDATED: v5, 2/21/23 - makes it so that it handles gff files insetad of feature tables because it appears that NCBI has depreciated these
UPDATED: v4, 12/13/22 - adds in annotation for the strand the hits fall on
UPDATE: v3, 2/11/21 - adds in annotation for "opposite to" genes.
This code is used to annotated the .tab results file of a cmsearch.
##########################
USEAGE:
MAKE SURE THAT YOUR FEATURE TABLES ARE SORTED WITH "removeduplicateentires.py" So that there are NO DUPLICATE ENTIRES.

How to prepare for running this code:

1) Run a cmsearch, and save your results as a .tab file. This is the "outputfromsearch" value.

2) Compile a large .txt file of all the gff3 from the NCBI FTP website. This can be done with the conda NCBI datasets download, where you can
select the gff3 files and concatentate them into a big .txt file.
Then run the "removeduplicateentires.py" script on this table. THIS IS IMPORTANT - If your genomes file contains duplicates the below code will not work properly!

3) Compile a large .fasta file of all the FASTA genomes (again, you can use ncbi datasets conda application for this). Concatenate these together.
Providing this FASTA file will annotate the results with the FASTA sequence for each hit. If you would like to skip this annotation step, write PASS instead of providing the FASTA file.

4) Name your output CSV file.
'''

import csv
import sys
import re
import pandas as pd
from Bio import SeqIO

##This is our main function
def annotatecmresults(outputfromsearch, gff3_reformatted, GenomesFASTA, outputfile):

##1) Read through the .tab output file and parse it so that all the spaces are "," between entries (I couldn't figure out a better way to do this!)"
    with open(outputfromsearch, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter='\t')
        genomenameslist =[]
        middlecoordinatelist = []
        strandlist = []
        startlist = []
        endlist = []

        for entry in fhcsv:
            entry = str(entry)
            entry = entry.rstrip("']")
            entry = entry.lstrip("['")
            entry = re.sub(' +', ',', entry)
            entry = entry.split(",")

            ##There are some random entries at the bottom of the .tab file but they start with #, so we are going to ignore those
            if "#" not in entry[0]:
                genomename = entry[0]
                start = entry[7]
                end = entry[8]
                middlecoordinate = round((int(start) + int(end)) / 2)
                strand = entry[9]

                ##Store all these important things in lists
                genomenameslist.append(genomename)
                middlecoordinatelist.append(middlecoordinate)
                strandlist.append(strand)
                startlist.append(start)
                endlist.append(end)

        if len(genomenameslist) != len(middlecoordinatelist) != len(strandlist):
            print("ERROR in code block 1: The length of the gene names, middle coordinate and strand lists are not equal!")
            exit()

        #Adding in the FASTA sequence of each hit to the record using BioPython
        FASTAseq = []
        description_lst = []


        if GenomesFASTA == "PASS":
            FASTAseq.append("No FASTA provided")

        for value in range(0, len(genomenameslist)):
            recordtracker = 0
            for seq_record in SeqIO.parse(GenomesFASTA, "fasta"):
                if genomenameslist[value] == seq_record.id:
                    recordtracker = recordtracker + 1
                    myseq = seq_record.seq
                    if startlist[value] != "N/A":
                        start = int(startlist[value])
                        end = int(endlist[value])
                    description = seq_record.description
                    description_lst.append(description)
                    coordinates = str(start) + ".." + str(end)

                    if start > end:
                        seqslice = myseq[end:start]
                        RCseqslice = seqslice.reverse_complement()
                        record = (">" + description + ", c" + coordinates + "\n" + RCseqslice)
                        FASTAseq.append(record)

                    if end > start:
                        seqslice = myseq[start:end]
                        record = (">" + description + ", " + coordinates + "\n" + seqslice)
                        FASTAseq.append(record)

            if recordtracker > 1:
                print("ERROR in code block 4: MORE THAN ONE GENOME WITH THAT NAME FOUND IN INPUT FASTA FILE.")

            if recordtracker == 0:
                record = (">ERROR: NO SEQ FOUND")
                FASTAseq.append(record)

            #Have to do this so that the fasta sequences appear correctly in the pandas df
            FASTAseq_updated = []
            for entry in FASTAseq:
                entry = str(entry)
                FASTAseq_updated.append(entry)

        #Lets turn this into a pandas dataframe
        data = [genomenameslist,middlecoordinatelist,strandlist, startlist, endlist, FASTAseq_updated, description_lst]
        cmresults_df = pd.DataFrame(data)
        cmresults_df = cmresults_df.transpose()
        cmresults_df.columns = (['Genome', "Mid_Coordinate", "Strand", "Start", "End", "FASTAseq", "Description"])


##2) Lets open the pre-formatted db as a pandas df
    genome_db_df = pd.read_csv(gff3_reformatted)


##3) Now we compare the lists - from the .tab results, and from the larger genome file.
    Annotated_cm_results_df = pd.DataFrame()

    cmresults_df_len = len(cmresults_df)
    for value in range (0, cmresults_df_len):

        #Select out the rows from the larger genome_db_df that only relate to the genome of the current cm hit we are iterrating through
        genome_name = cmresults_df.iloc[value, 0] #This is the genome name of the current cm result
        hit_mid_cord = cmresults_df.iloc[value, 1]
        hit_strand = cmresults_df.iloc[value, 2]
        hit_start = cmresults_df.iloc[value, 3]
        hit_end = cmresults_df.iloc[value, 4]
        hit_fasta = cmresults_df.iloc[value, 5]
        genome_description = cmresults_df.iloc[value, 6]

        selected_genome_db_df = genome_db_df[genome_db_df["Genome"]== genome_name] #Selecting out the rows from the bigger db file that have the same names as the
        selected_genome_db_df_len = len(selected_genome_db_df)

        for row in range(0, selected_genome_db_df_len):
            Start = selected_genome_db_df.iloc[row, 4]
            End = selected_genome_db_df.iloc[row, 5]

            if hit_mid_cord in range (Start, End):
                #If the result is in an annotated gene, we want to store all the info about that hit, and the genes left/right of it.
                Genome = selected_genome_db_df.iloc[row, 0]
                Gene = selected_genome_db_df.iloc[row, 1]
                FeatureType = selected_genome_db_df.iloc[row, 2]
                Strand = selected_genome_db_df.iloc[row, 3]
                Locus_tag = selected_genome_db_df.iloc[row, 6]
                Product = selected_genome_db_df.iloc[row, 7]

                #I want to note if the hit and the gene are on the same strand, or opposite
                if str(FeatureType) == "IGR":
                    Strand_match = str("IGR")

                if str(Strand) == "-" and str(hit_strand) == "plus":
                    Strand_match = str("Opposite")

                if str(Strand) == "+" and str(hit_strand) == "minus":
                    Strand_match = str("Opposite")

                if str(Strand) == "+" and str(hit_strand) == "plus":
                    Strand_match = str("Same")

                if str(Strand) == "-" and str(hit_strand) == "minus":
                    Strand_match = str("Same")

                #Gene left
                left_gene = int(row -1) #I don't think we need to worry about the first entry (and last for the right gene), becuase the genome is circular. So the -1 entry is the last entry, which is fine.

                #if the gene left is an IGR, we don't want that info! We only want actual annotations!
                left_gene_FeatureType = selected_genome_db_df.iloc[left_gene, 2]
                if left_gene_FeatureType =="IGR":
                    left_gene = int(row -2)

                #Save all the info from the left gene
                left_gene_Gene = selected_genome_db_df.iloc[left_gene, 1]
                left_gene_FeatureType = selected_genome_db_df.iloc[left_gene, 2]
                left_gene_Strand =selected_genome_db_df.iloc[left_gene, 3]
                left_gene_Start = selected_genome_db_df.iloc[left_gene, 4]
                left_gene_End = selected_genome_db_df.iloc[left_gene, 5]
                left_gene_Locus_tag = selected_genome_db_df.iloc[left_gene, 6]
                left_gene_Product = selected_genome_db_df.iloc[left_gene, 7]

                #Gene right
                right_gene = int(row +1)

                #if the gene right is an IGR, we don't want that info! We only want actual annotations!
                right_gene_FeatureType = selected_genome_db_df.iloc[right_gene, 2]
                if right_gene_FeatureType =="IGR":
                    right_gene = int(row + 2)

                #Save all the info from the right gene
                if right_gene < selected_genome_db_df_len: #if the hit is the last gene in the genome, the gene right is actually going to be the first gene in the genome (because it's circular)
                    right_gene_Gene = selected_genome_db_df.iloc[right_gene, 1]
                    right_gene_FeatureType = selected_genome_db_df.iloc[right_gene, 2]
                    right_gene_Strand =selected_genome_db_df.iloc[right_gene, 3]
                    right_gene_Start = selected_genome_db_df.iloc[right_gene, 4]
                    right_gene_End = selected_genome_db_df.iloc[right_gene, 5]
                    right_gene_Locus_tag = selected_genome_db_df.iloc[right_gene, 6]
                    right_gene_Product = selected_genome_db_df.iloc[right_gene, 7]

                if right_gene > selected_genome_db_df_len: #This is in that special case where the hit was at the end of the genome, and the right gene would then be the first gene in the genome
                #I don't know if this will handle IGRs, but I'm too lazy right now to accomidate this super specific and special case. If it happens to you, sorry. But the start of the genome probably shouldn't be an IGR...
                    right_gene_Gene = selected_genome_db_df.iloc[0, 1]
                    right_gene_FeatureType = selected_genome_db_df.iloc[0, 2]
                    right_gene_Strand =selected_genome_db_df.iloc[0, 3]
                    right_gene_Start = selected_genome_db_df.iloc[0, 4]
                    right_gene_End = selected_genome_db_df.iloc[0, 5]
                    right_gene_Locus_tag = selected_genome_db_df.iloc[0, 6]
                    right_gene_Product = selected_genome_db_df.iloc[0, 7]

                #This new_row variable is where we're going to put all the info that we want to annotate into the final results!
                new_row = pd.DataFrame({
                "Genome":[Genome],"Description":[genome_description],

                "left_gene_Gene":[left_gene_Gene], "left_gene_FeatureType":[left_gene_FeatureType], "left_gene_Strand":[left_gene_Strand],
                "left_gene_Start":[left_gene_Start], "left_gene_End":[left_gene_End], "left_gene_Locus_tag":[left_gene_Locus_tag], "left_gene_Product":[left_gene_Product],

                "Annotated_Gene":[Gene], "Annotated_FeatureType":[FeatureType], "Annotated_Strand":[Strand],
                "Annotated_Start":[Start],"Annotated_End":[End], "Annotated_Locus_tag":[Locus_tag], "Annotated_Product":[Product],
                "cm_Hit_Strand":[hit_strand], "Orientation_to_annotation":[Strand_match], "cm_Hit_Start":[hit_start], "cm_Hit_End":[hit_end], "cm_Hit_FASTA":[hit_fasta],

                "right_gene_Gene":[right_gene_Gene], "right_gene_FeatureType":[right_gene_FeatureType], "right_gene_Strand":[right_gene_Strand],
                "right_gene_Start":[right_gene_Start], "right_gene_End":[right_gene_End], "right_gene_Locus_tag":[right_gene_Locus_tag], "right_gene_Product":[right_gene_Product]
                })

                Annotated_cm_results_df = pd.concat([Annotated_cm_results_df, new_row], join="outer")

    #Saving it to a file
    Annotated_cm_results_df.to_csv(outputfile, index=False)


if __name__ == '__main__':
    if len(sys.argv) == 5:
        annotatecmresults(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("1) the .tab output from the cmsearch; 2) annoated gff3 file as a csv; 3) a large FASTA file of all the genomes - if you want to pass this option write PASS for this value; 4) Output file name (CSV)")
        sys.exit(0)
