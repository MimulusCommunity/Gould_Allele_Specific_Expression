## This script will take the output of two blast searches and return the reciprocal best (first) hits in list form with alignment ranges from the best alignment.

## This script is designed to work where the second input file is from a transcriptome with multiple alleles/pair. RBB hits are returned at the level of the pair/gene.  ##

# input file format:
# 0:qseqid 1:sseqid 2:pident 3:length 4:mismatch 5:gapopen 6:qstart 7:qend 8:sstart 9:send 10:evalue 11:bitscore 12:qcovhsp
Usage = "Useage: script Ref_vs_transcriptome.out Transcriptome_vs_ref.out > RBH-list-outfile &2>run.log "

import sys, re

if len(sys.argv) < 3:
	print(Usage)

def blast_list_parser(open_file):
	D = {} #dictionary for BLAST file
	# prev_gene_data = [queryID, queryId_alnA_start, queryId_alnA_end, subjectId, subjectId_alnA_start, subjectId_alnA_end, alnA_eval, pct_cov, num_mismatches, num_gaps]
	hsp_list = [ ("#", 0, 0, "", 0, 0, 1, 0) ]
	num_mult_hsps = 0
	for Line in open_file:
		if ( "#" not in Line):
		#	print(Line)
			data=Line.strip()
			Elements = re.split('\t', data)
			line_data = (re.split("_",Elements[0][0:12])[0], Elements[6], Elements[7], re.split("_", Elements[1][0:12])[0], Elements[8], Elements[9], Elements[10], Elements[12], Elements[4],Elements[5])
			queryId = re.split("_", Elements[0][0:12])[0]	
		#	queryId_alnA_start = Elements[6]
		#    queryId_alnA_end = Elements[7]
		#	subjectId = Elements[1]
		#    subjectId_alnA_start = Elements[8]
		#    subjectId_alnA_end = Elements[9]
			#sys.stderr.write( "QID: " + str(queryId) + " LastHSP: " +str(hsp_list[-1][0]) )
			if queryId == hsp_list[-1][0]:			# if the present gene matches the previous gene (is a multiple HSP), add dataline to list
				hsp_list.append(line_data)

			if queryId != hsp_list[-1][0]:	# if the present gene doesnt match the last hsp
				
				if len(hsp_list) > 1: num_mult_hsps += 1
				sum_hsp_cov = float(0.0) 			# sum the coverage for the prev hsps
				perc_cov_list = []
				
				for dataline in hsp_list:
					perc_cov = int(dataline[7])
					perc_cov_list.append(perc_cov)
					sum_hsp_cov = sum_hsp_cov + perc_cov
				if sum_hsp_cov >= 75:			#if the sum is greater than 75%, add longest coverage hsp to dict
					index = perc_cov_list.index(max(perc_cov_list))
					D[hsp_list[index][0]] = hsp_list[index][1:]
				hsp_list = [(line_data)]		# start a new hsp list with the current data line
		
	if len(hsp_list) > 0:					#write the last HSP
		sum_hsp_cov = int(0)                         # sum the coverage for the prev hsps
                perc_cov_list = []
                for dataline in hsp_list:
                       perc_cov = int(dataline[7])
                       perc_cov_list.append(perc_cov)
                       sum_hsp_cov = sum_hsp_cov + perc_cov
                       if sum_hsp_cov >= 75:                   #if the sum is greater than 75%, add longest coverage hsp to dict
                                index = perc_cov_list.index(max(perc_cov_list))
                                D[hsp_list[index][0]] = hsp_list[index][1:]

	sys.stderr.write("Parsed a file.\n")
	return (D, num_mult_hsps)

#########################################
sys.stderr.write(" \n \n Input files must be sorted by queryID.\n\n")

infl1 = sys.argv[1]
infl2 = sys.argv[2]

#parse first BLAST results
fileA = open(infl1, 'r')
(D1,N) = blast_list_parser(fileA)
sys.stderr.write("Num_mult_hsps: " + str(N) +"\n")
fileA.close()

#sys.exit()

#parse second BLAST results with "B" alignments
fileB = open(infl2, 'r')
(D2,N) = blast_list_parser(fileB)
sys.stderr.write("Num_mult_hsps: " + str(N) +"\n")
fileB.close()

#Now, pick the share pairs
# D[gene] = [ queryId_alnA_start, queryId_alnA_end, subjectId, subjectId_alnA_start, subjectId_alnA_end, alnA_eval, pct_cov]

sys.stderr.write("D1: " + str(D1)[0:100])
sys.stderr.write("D2: " + str(D2)[0:100])

SharedPairs={}
for gene_name in D1.keys():
#	print("testing: " + gene_name)
	match_name = D1[gene_name][2].split("_")[0]     #gene1 is matched by GENE_PAIR e.g. "pair10_"
#	print("match name: " + match_name)

	matches = []
	evals = []
	for item in D2.keys():							# find multiple alleles for the GENE_PAIR
		gene_pair = item.split("_")[0]
		#print(gene_pair)
		if (gene_pair == match_name) & (str(D2[item][2]) == gene_name):					# this equality is a RBB Hit.
			matches.append(item)
			evals.append(D2[item][5])
		
#			print("Matches: " + str(matches) + "\n")
#			print("Evals: " + str(evals) + "\n")

	if len(matches) >= 1:
		best_match = matches[evals.index(min(evals))]			# find the best matching allele via evalue
		SharedPairs[gene_name] = ('\t'.join(D1[gene_name][0:2]), best_match, '\t'.join(D2[best_match][0:2]))		# capture the alignment scores for the best matching allele
#		print("Best match: " + str(best_match) + "\t"+ str(SharedPairs[gene_name]))
#	else:
#		print("No match")

sys.stdout.write("RefGene\tg1_start\tg1_end\tAlleleMatch\tg2_start\tg2_end\n")

count = 0
for entry in SharedPairs.keys(): 
	line = str("\t".join(SharedPairs[entry]))
	sys.stdout.write(entry + "\t" + line + "\n")
	count += 1

sys.stderr.write("There are: " + str(len(SharedPairs.keys())) + " RBB pairs.\n")
#sys.stderr.write("RBB Pairs: " + str(count) + "\n")
sys.stderr.write("COMPLETE")
