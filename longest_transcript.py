# This script will take a list of contig alignments and output allele pair number and alignment length.

import sys

alignments = open(sys.argv[1], "r")

pair_no = 0

for line in alignments:
    if "gene1" not in line:
        pair_no += 1
        start1, end1, start2, end2 = map(int, [ line.split()[i] for i in [1,2,4,5]])
        len1 = max(start1, end1) - min(start1, end1)
        len2 = max(start2, end2) - min(start2, end2)
        sys.stdout.write("pair" + str(pair_no) + "\t" + str(max(len1, len2)) + "\n")

alignments.close()

sys.stderr.write("Complete")

