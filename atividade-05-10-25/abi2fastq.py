from Bio import SeqIO
import sys

input = sys.argv[1]

out = sys.argv[2]

records = SeqIO.parse(input, "abi")
count = SeqIO.write(records, out, "fastq")
print("Converted %i records" % count)
