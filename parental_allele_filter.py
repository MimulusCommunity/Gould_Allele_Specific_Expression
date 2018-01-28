# This script will rewrite parental genotype data for all genes with allele ratio above a certain threshold
# *** input file must be sorted by gene and then by allele

import sys, re

####################
# a function for outputing the allele ratios for a parental allele count file. examine histogram
#infile = open(sys.argv[1], "r")

#prev_gene="x"

#sys.stdout.write("Gene \t L1/S1 ratio \n")

#for line in infile:
#    if "Gene" not in line:
#        gene = line.split()[0]
#
#        if gene != prev_gene:
#            Lallele = float(line.split()[5])
#            prev_gene = gene
#
#        elif gene == prev_gene:
#            Sallele = float(line.split()[5])
#            if Sallele == 0:
#                ratio = 0
#            else:
#                ratio = float(Lallele/Sallele)
#            sys.stdout.write(gene + "\t" + str(ratio) + "\n")

#infile.close()
#print("Complete")
#######################



#############
# a function for identifying all genes with low signal to noise ratio

# input *.counts file. FORMAT: Gene \t Allele1 \t Allele2

library = str(sys.argv[1])
params = re.split(r'_|-', library)    # (25, P6, 3 , 24, L1, 12, CGCTCATT, L006.counts)
parent = params[4]
elim = 0

#if parent != "L1" and parent != "S1":
#    sys.stderr.write("Err: This is not a parental count file.")

#infile = open(sys.argv[1], "r")
#for line in infile:
    #print(line)
#    if "pair" in line:              # skip header rows
#        gene = line.split()[0]
#	print(str(gene))

#        if parent == "L1":
#            parent_ac = float(line.split()[1])
#            alt_ac = float(line.split()[2])
#            print(str(parent_ac))
       
#	elif parent == "S1":
#	    parent_ac = float(line.split()[2])
#            alt_ac = float(line.split()[1])
#            print(str(alt_ac))
  # if the correct/incorrect allele ratio is <5,  print to elim list
 #       if float(parent_ac/5) < alt_ac:
 #           elim += 1
#	    sys.stdout.write(gene +"\n")

#infile.close()
#sys.stderr.write("Library: " + library + "  Eliminated genes: " + str(elim) + "\n")



############
# a function to filter out low signal to noise ratio genes

genes = open(sys.argv[2], "r")
elim_list=[]
for gene in genes:
    elim_list.append(gene.strip())
    #sys.stderr.write(str(elim_list[0:10]))
genes.close()

infile = open(sys.argv[1], "r")
for line in infile:
    if "pair" in line:              # skip header
        gene = line.split()[0]
	#sys.stderr.write(gene)
        if gene not in elim_list:
	    sys.stdout.write(line)
infile.close()
sys.stderr.write(library + " complete. \n")
