from Bio import SeqIO
import sys

input_file = sys.argv[1]
cutoff_length = int(sys.argv[2])

fasta_sequences = SeqIO.parse(open(input_file),'fasta')

for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    seqlen = len(sequence)
    new_sequence = str(sequence).split("*")
    new_sequence = new_sequence[0]
    if(name.endswith("_1")):
        name = name.replace("|","_1|")
    if(name.endswith("_2")):
        name = name.replace("|","_2|")
    if(name.endswith("_3")):
        name = name.replace("|","_3|")
    if(name.endswith("_4")):
        name = name.replace("|","_4|")
    if(name.endswith("_5")):
        name = name.replace("|","_5|")
    if(name.endswith("_6")):
        name = name.replace("|","_6|")
    if (len(new_sequence)>cutoff_length):
        print(">"+str(name)+"\n"+new_sequence)
