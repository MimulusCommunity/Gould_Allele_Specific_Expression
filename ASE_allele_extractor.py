# This script will take a list of transcripts from two parental transcriptomes and extract them at the aligned regions into a combined reference fasta.

import sys

def parse_transcriptome(infile, numseqs, prefix):
    D = {}
    for i in range(0, numseqs):
        line = infile.readline().split()
	if ">" in line[0]:
	    contig = (str(prefix) + str(line[0][1:]))
            sequence = infile.readline()
            D[contig] = sequence
    sys.stderr.write("Parsed a file: " + str(len(D.keys())) + " contigs.\n")
    #print(D.keys()[0:5])
    return(D)

###################

# read tr1 into a dictionary.
transcriptome1 = open(sys.argv[1], "r")
numseqs = len(transcriptome1.read().splitlines())/2
transcriptome1.close()

transcriptome1 = open(sys.argv[1], "r")
DL1 = parse_transcriptome(transcriptome1, numseqs, "L1-")
transcriptome1.close()


#read tr2 into a dictionary.
transcriptome2 = open(sys.argv[2], "r")
numseqs = len(transcriptome2.read().splitlines())/2
transcriptome2.close()

transcriptome2 = open(sys.argv[2], "r")
DS1 = parse_transcriptome(transcriptome2, numseqs, "S1-")
transcriptome2.close()

#for each line in the Recip Best Blast list, extract the aligned regions, look up sequences, extract aligned regions into a file.

# FORMAT:
# gene1	g1_start	g1_end	gene2	g2_start	g2_end	e-val	pct_cov
# L1-c67197_g16_i1	39	566	S1-c39672_g2_i7166	694	0.0	88
# S1-c62991_g1_i1	1	590	L1-c78413_g1_i1	39	628	0.0	100

def extract_contig_reg(contig, start, end):
    #print(contig)
    real_start = min(start,end)     # for rev. comp. alignments, the downstream base is listed first in the alignment file, so this has to be taken into account.
    real_end = max(start,end)
    if contig.startswith("L1-"):
        sequence = DL1[contig]
        region = sequence[(real_start-1):real_end]
    elif contig.startswith("S1-"):
        sequence = DS1[contig]
    	region = sequence[(real_start-1):real_end]

    return(region)

pairs = open(sys.argv[3], "r")
pair_number = 0
for line in pairs:
    if "gene1" not in line:
        pair_number += 1
    
        contig1 = line.split()[0]
        contig1_start = int(line.split()[1])
        contig1_end = int(line.split()[2])
        region1 = extract_contig_reg(contig1, contig1_start, contig1_end)
        contig1_name = "pair" + str(pair_number) + "_" + contig1
        sys.stdout.write( ">" + str(contig1_name) + "\n" + str(region1) + "\n")
    
        contig2 = line.split()[3]
        contig2_start = int(line.split()[4])
        contig2_end = int(line.split()[5])
        region2 = extract_contig_reg(contig2, contig2_start, contig2_end)
        contig2_name = "pair" + str(pair_number) + "_" + contig2
        sys.stdout.write( ">" + str(contig2_name) + "\n" + str(region2) + "\n")

sys.stderr.write("Complete. Pairs written to combined reference: " + str(pair_number) + "\n")

