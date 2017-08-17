from Bio import SeqIO
import sys

input_file = sys.argv[1]
length = int(sys.argv[2])

fasta_sequences = SeqIO.parse(open(input_file),'fasta')

counter = 0
for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    if(len(sequence)<=length):
        counter += 1
print(counter)
