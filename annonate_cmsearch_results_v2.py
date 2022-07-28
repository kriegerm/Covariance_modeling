#!/usr/bin/env python3

'''
annonate_cmsearch_results_v2.py
UPDATE: v2, 1/12/20
This code is used to annotated the .tab results file of a cmsearch.

##########################
USEAGE:

MAKE SURE THAT YOUR FEATURE TABLES ARE SORTED WITH "removeduplicateentires.py" So that there are NO DUPLICATE ENTIRES (see note 2 below). The output file must be in the following order:

How to prepare for running this code:
1) Run a cmsearch, and save your results as a .tab file. This is the "outputfromsearch" value.

2) Compile a large .txt file of all the Feature Tables from the NCBI FTP website. This can be done by running a shell script with the FTP addresses to these feature tables, as below:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/267/845/GCF_001267845.1_ASM126784v1/GCF_001267845.1_ASM126784v1_feature_table.txt.gz
One you download and unzip all these .txt feature tables, concatenate them into one big file.
Then run the "removeduplicateentires.py" script on this table. THIS IS IMPORTANT - If your tabular genomes file contains duplicates the below code will not work properly!

3) Compile a large .fasta file of all the FASTA genomes from the NCBI FTP website. You can run a similar shell script to get these, and then concatenate all together:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/267/845/GCF_001267845.1_ASM126784v1/GCF_001267845.1_ASM126784v1_genomic.fna.gz
Providing this FASTA file will annotate the results with the FASTA sequence for each hit. If you would like to skip this annotation step, write PASS instead of providing the FASTA file.

4) Name your output CSV file.

'''

import csv
import sys
import re
from Bio import SeqIO


##This is the function we will use to make the list of indices in code block 3
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


##This is our main function
def annotatecmresults(outputfromsearch, FeatureTables, GenomesFASTA, outputfile):

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


##2) Following the same approach, we make a similar lists of all the items in the concatenated tabular genomes file that you downloaded from the fth website on NCBI
##THIS FILE HAS TO BE EDITED WITH THE "removeduplicateentries.py" CODE!
    with open(FeatureTables, 'r') as tabfh:
        tabfhcsv = csv.reader(tabfh, delimiter=",")
        field_names_list = next(tabfhcsv)

        db_genomename = []
        db_featuretype = []
        db_genestrand = []
        db_genename = []
        db_genestart =[]
        db_geneend = []
        db_genesymbol = []

        for line in tabfhcsv:
            if line[0] != "Genome":
                featuretype = line[3]  # gene, CDS, tRNA, etc
                genomeassembly = line[0]  # name of organism, need to match this to the genomename from the .tab output
                genestart = line[4]
                geneend = line[5]
                genestrand = line[6]
                genename = str(line[1])
                genesymbol = str(line[2])

                db_genename.append(genename)
                db_genomename.append(genomeassembly)
                db_featuretype.append(featuretype)
                db_genestrand.append(genestrand)
                db_genestart.append(genestart)
                db_geneend.append(geneend)
                db_genesymbol.append(genesymbol)

        if len(db_genename) != len(db_geneend) != len(db_genestart) != len(db_genestrand) != len(db_featuretype) != len(db_genomename) != len(db_genesymbol):
            print("ERROR in code block 2: the length of the gene name list, gene start and end lists, gene strand, genome name, and feature type lists are not equal!")
            exit()


##3) Now we compare the lists - from the .tab results, and from the larger genome file.
        gene = []
        symbol = []
        gene_start = []
        gene_end = []
        genome = []
        left = []
        right = []


        for value in range(0, len(genomenameslist)): #this is correct

            intmiddle = int(middlecoordinatelist[value])
            indices = get_index_positions(db_genomename, genomenameslist[value])            # this returns a list of all the indices of where the specific result is occuring in the larger tabular list

            if not indices:
                gene.append("N/A")
                genome.append("No genome found for this search")
                left.append("N/A")
                right.append("N/A")
                gene_start.append("N/A")
                gene_end.append("N/A")
                symbol.append("N/A")
                continue

            count = 0

            for location in indices:

                intstart = int(db_genestart[location])
                intend = int(db_geneend[location])

                if intmiddle in range(intstart, intend):
                    count = count +1

                    gene.append(db_genename[location])
                    genome.append(db_genomename[location])
                    gene_start.append(db_genestart[location])
                    gene_end.append(db_geneend[location])
                    symbol.append(db_genesymbol[location])

                    L1 = int(int(location) -1)
                    L2 = int(int(location) -2)
                    R1= int(int(location) + 1)
                    R2 = int(int(location) + 2)


                    if str(db_genename[L1]) == "IGR":
                        left.append(db_genename[L2])

                    if str(db_genename[L1]) != "IGR":
                        left.append(db_genename[L1])

                    if str(db_genename[R1]) == "IGR":
                        right.append(db_genename[R2])

                    if str(db_genename[R1]) != "IGR":
                        right.append(db_genename[R1])

                    continue


            if count == 0: #if there are no hits where the middle value of the cooridnates falls between two genome locations, then the following values will be appended to the list:
                genome.append(genomenameslist[value])
                gene.append(genomenameslist[value])
                left.append("None found - coordinates out of range")
                right.append("None found - coordinates out of range")
                gene_start.append("None found - coordinates out of range")
                gene_end.append("None found - coordinates out of range")
                symbol.append("None found - coordinates out of range")

            if count > 1:
                print("ERROR in code block 3: more than one entry found for a result in the .tab file!")
                print("The offending entry is: " + str(genomenameslist[value]))
                exit()


        if len(gene) != len(genome) != len(right) != len(left):
            print("ERROR in code block 3: The number of genes and left/right matches are not equal!")
            exit()

        if len(genome) != len(genomenameslist):
            print("ERROR in code block 3: different numbers of genes found than results in the .tab file.")
            print(genome)
            print(genomenameslist)
            exit()


##4) Add the FASTA Seqnece of the hit to the file - if the input value for this is PASS then it is not annotated
        FASTAseq = []

        if GenomesFASTA == "PASS":
            FASTAseq.append("No FASTA provided")

        for value in range(0, len(genomenameslist)):
            recordtracker = 0

            for seq_record in SeqIO.parse(GenomesFASTA, "fasta"):


                if genome[value] == seq_record.id:
                    recordtracker = recordtracker + 1
                    myseq = seq_record.seq

                    if gene_start[value] != "N/A":
                        start = int(startlist[value])
                        end = int(endlist[value])

                    description = seq_record.description
                    coordinates = str(start) + ".." + str(end)

                    if start > end:
                        # if strand == "minus":
                        seqslice = myseq[end:start]
                        RCseqslice = seqslice.reverse_complement()
                        record = (">" + description + ", c" + coordinates + "\n" + RCseqslice)
                        FASTAseq.append(record)

                    if end > start:
                         # if strand == "plus":
                        seqslice = myseq[start:end]
                        record = (">" + description + ", " + coordinates + "\n" + seqslice)
                        FASTAseq.append(record)

            if recordtracker > 1:
                print("ERROR in code block 4: MORE THAN ONE GENOME WITH THAT NAME FOUND IN INPUT FASTA FILE.")
                exit()

            if recordtracker == 0:
                record = (">ERROR: NO SEQ FOUND")
                FASTAseq.append(record)


##5) Write the results file
    with open(outputfile, 'a') as fh:
        writer = csv.writer(fh)
        header = ["Genome", "Gene Left", "Hit Annotation", "Hit Symbol", "Gene Right", "FASTA Sequence of Hit"]
        writer.writerow(header)
        for value in range(0, len(genomenameslist)):
            saver = []
            saver.append(genome[value])
            saver.append(left[value])
            saver.append(gene[value])
            saver.append(symbol[value])
            saver.append(right[value])
            saver.append(FASTAseq[value])
            writer.writerow(saver)


if __name__ == '__main__':
    if len(sys.argv) == 5:
        annotatecmresults(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("1) the .tab output from the cmsearch; 2) a large file of all the feature_table.txt that is found on the NCBI FTP website with NO DUPLICATES AND IGRS ANNOTATED; 3) a large FASTA file of all the genomes from the NCBI FTP website - if you want to pass this option write PASS for this value; 4) Output file name (CSV)")
        sys.exit(0)

