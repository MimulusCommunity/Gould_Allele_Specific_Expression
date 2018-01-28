# This script is a simple counting script for reads aligned to a transcriptome. Raw read counts for each contig are returned.
# Data output is formatted for allele specific expression analysis.

## Useage: script in.sam ref.fa > out.txt

import sys, re

# Sample naming convention:
#L1-	inland
#S1 	coastal
#SL	F1
#0-100	bodega
#101-200	pepperwood
#P	pool number

# set the rep#, environment, and gneration values according to the library name
# ex. 25_P6-3_24-SL-12_CGCTCATT_L006.sam

library = str(sys.argv[1])
params = re.split(r'_|-', library)    # (25, P6, 3 , 24, SL, 12, CGCTCATT, L006.sam)

Rep = params[5]

if params[4] == "SL":   # Set the Generation
    Gen = "F1"
elif params[4] == "S1" or params[4] == "L1":
    Gen = "P"
else:
    sys.stderr.write("Unknown Gen: " + str(params[4]) +"\n")

if int(params[5]) in range(0,100):      # Set the Environment
    Env = "C"
elif int(params[5]) in range(101,200):
    Env = "I"
else:
    sys.stderr.write("Unknown Env: " + str(params[5] + "\n"))



# make a dictionary of all gene/pair names from the reference file:

ref = open(sys.argv[2], "r")

Counts = {}

for line in ref:
    if ">" in line:
        name = line.split()[0][1:]
#	allele = re.split(r'-',name)[0]
        pair = re.split(r'_|-', name)[0]        # pair14413
#        parent_allele = re.split(r'_|-', name)[1]      # L1
        Counts[pair] = [pair, str(0), str(0)]      # set count L1 = 0, S1 = 0

sys.stderr.write("There are " + str(len(Counts.keys())) + " contigs in the reference \n Counting reads . . .\n")

#create dictionary of counts per gene output in following format
# 
# Example: pair14413_L1-c66299_g1_i12

sam_file = open(sys.argv[1], "r")

unaligned = 0
aligned = 0

for line in sam_file:
    if "@" not in line:
	try:
            contig = line.split()[2]        # pair14413_L1-c66299_g1_i12
	except IndexError:
	    sys.stderr.write("Bad line: " + line)
	    continue				#this happend very very occasionally due to errors in the SAM files, I think
 
        if contig == "*":
            unaligned += 1

        else:
	    aligned += 1
#	    name = re.split(r'-', contig)[0]
	    try:
	        pair = re.split(r'-|_', contig)[0]
                parental_allele = re.split(r'_|-', contig)[1]
	    except IndexError:
		sys.stderr.write("Bad line: " + line)
		continue

	    if parental_allele == "L1":
            	prev_count = int(Counts[pair][1])
		count = prev_count + 1
                Counts[pair][1] = str(count) 
            elif parental_allele == "S1":
		prev_count = int(Counts[pair][2])
	        count = prev_count + 1
                Counts[pair][2] = str(count)

sam_file.close()

sys.stdout.write("Allele\t" + "L" + "\t" + "S" + "\n")
sys.stdout.write("Gen\t" + str(Gen) + "\t" + str(Gen) + "\n")
sys.stdout.write("Env\t" + str(Env) + "\t" + str(Env) + "\n")

for item in Counts.keys():
    sys.stdout.write('\t'.join(Counts[item]) + "\n")

sys.stderr.write("Count complete. \n Aligned reads: " + str(aligned) + "\n" + "Unaligned reads: " + str(unaligned) + "\n")

