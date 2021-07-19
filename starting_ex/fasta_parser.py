from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import transcribe 
from Bio.Seq import Seq 

file = open("human_mx1.fasta") 

#get the first itterable with next() 
sequence_record = next(parse(file, 'fasta'))

#print(sequence_record.id)

rna_sequence = transcribe(sequence_record.seq)
rna_sequence.back_transcribe()

protein_sequence = rna_sequence.translate()

print(protein_sequence)





















#old bits
""" records = parse(file, 'fasta') 
for record in records:
    print("Id: %s" % record.id) 
    print("Name: %s" % record.name) 
    print("Description: %s" % record.description) 
    print("Annotations: %s" % record.annotations) 
    print("Sequence Data: %s" % record.seq)
    dna_seq = record.seq

#print(dna_seq)

rna_seq = transcribe(dna_seq)
rna_seq.back_transcribe()

#print(rna_seq)

protein_seq = rna_seq.translate()

print(protein_seq)
 """
