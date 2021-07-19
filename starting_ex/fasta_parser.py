from Bio.SeqIO import parse, to_dict 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import transcribe 
from Bio.Seq import Seq 
from Bio.Blast import NCBIWWW

file = open("human_mx1.fasta") 

#get the first itterable with next() 
#sequence_record = next(parse(file, 'fasta'))

sequence_record_dict = to_dict(parse(file, 'fasta'))

#convert FIRST (0) dictionary keys to a list 
sequence_record = sequence_record_dict[list(sequence_record_dict)[0]]





#print(sequence_record.id)

rna_sequence = transcribe(sequence_record.seq)
rna_sequence.back_transcribe()

protein_sequence = rna_sequence.translate()

print(protein_sequence)

#result_handle = NCBIWWW.qblast("blastn", "nt", sequence_record.seq)



















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
