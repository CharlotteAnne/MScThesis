from Bio import SeqIO
import sys

input_file = sys.argv[1]

fasta_sequences = SeqIO.parse(open(input_file),'fasta')

for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    name = name.replace("|","_retained_intron|")
    print(">"+str(name)+"\n"+sequence)
