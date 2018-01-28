# This script will take a list of allele best blast hits, output the gene, hit annotation, inversion status, and reference coordinates in the genome.

# Useage : script.py blastoutput.txt CR.fa ref.gff > annotation.txt  *** Input files must be sorted ***

import sys
import re

# read in gff coordinate annotations as dictionary and annotate inversion genes

# GFF file format

###gff-version 3
##annot-version v2.0
#scaffold_1	phytozomev10	gene	16375	23539	.	-	.	ID=Migut.A00001.v2.0;Name=Migut.A00001;ancestorIdentifier=mgv1a001726m.g.IM62.v1.1
#scaffold_1	phytozomev10	mRNA	16375	23539	.	-	.	ID=Migut.A00001.1.v2.0;Name=Migut.A00001.1;pacid=28939341;longest=1;Parent=Migut.A00001.v2.0

RBBHs = open(sys.argv[1], "r")
D_RBBHs={}
for line in RBBHs:
    gene = line.split()[0]
    pairmatch = line.split()[3]
    D_RBBHs[gene] = pairmatch
RBBHs.close()

gff = open(sys.argv[2], "r")
D_gff = {}

for line in gff:
    if "#" not in line and "gene" in line:
        name = (line.split()[8])[3:15]
        chr = line.split()[0]
        start = (line.split()[3])
        end = (line.split()[4])
        strand = line.split()[6]

        if chr == "scaffold_5":              #inv 5 scaffolds
            left=(11466582,11781205,12467777,13082240,13367321,14418230,15044886,16114744,16561107,24038195)
            right=(11646589,12457777,12605828,13357321,13932160,14901369,15663194,16283041,18317358,24294111)
            for i in range(0,10):
                if int(left[i]) < int(start) < int(right[i]):
                    inversion = "I5"

        elif chr == "scaffold_8" and (888603 < int(start) < 7594555):    # scaffold 8 from 888603 to 7594555
            inversion = "I8"

        else:
            inversion = "NON"
        
        if name in D_RBBHs.keys():
            pairmatch = D_RBBHs[name]
        else:
            pairmatch = "--"
            
        D_gff[name]=[pairmatch, chr, start, end, strand, inversion]
gff.close()

# blast output format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp

for gene in D_gff.keys():
    sys.stdout.write(gene +"\t"+ "\t".join(D_gff[gene]) +"\n")

sys.stderr.write( "Annotation complete.")

sys.stderr.write( "REMEMBER TO SORT THE OUTPUT FILE BY TRANSCRIPT \n")


