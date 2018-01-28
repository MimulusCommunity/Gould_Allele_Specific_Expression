## This script will filter out a list of sequences (contaminants) from a large fasta file, returning the cleaned file.
# Writes to std.out

# First remove the extra line breaks from your input sequence file! Use this AWK program:
# % awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' input.fasta > output.fasta

import sys

dirty_File = open(sys.argv[1], "r")
genelist=dirty_File.read().splitlines()
numlines = len(genelist)
dirty_File.close()
sys.stderr.write( "There are " + str(numlines/2) + " sequences in the original file \n")


contaminants = open(sys.argv[2], "r")
names_to_remove = []
for line in contaminants:
    name = line.split()[0]
    names_to_remove.append(name)
contaminants.close()

#sys.stderr.write(str(names_to_remove[0:8]))

sys.stderr.write("Removing " + str(len(names_to_remove)) + " sequences (including duplicates) . . .\n")

dirty_File = open(sys.argv[1], "r")
good_contigs = 0
for i in range(0, numlines):
    dataline = str( dirty_File.readline() )
    if ">" in dataline:
        contig_name = dataline.split()[0][1:]
#        sys.stderr.write(str(contig_name))
	if contig_name not in names_to_remove:
            sequence = dirty_File.readline()
            sys.stdout.write(">"+str(contig_name) + "\n" + str(sequence))
            good_contigs += 1

sys.stderr.write("Returning " + str(good_contigs) + " sequences.\nComplete.")

dirty_File.close()
