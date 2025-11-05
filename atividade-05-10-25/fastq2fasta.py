from Bio import SeqIO
import sys

input = sys.argv[1]

out = sys.argv[2]

records = SeqIO.parse(input, "fastq")
count = SeqIO.write(records, out, "fasta")
print("Converted %i records" % count)
