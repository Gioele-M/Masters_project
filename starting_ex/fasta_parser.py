import pandas as pd
import numpy as np
from Bio.SeqIO import parse, to_dict 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import transcribe 
from Bio.Seq import Seq 
from Bio.Blast import NCBIWWW
from Bio import SearchIO


file = open("human_mx1.fas") 

#get the first itterable with next() 
#sequence_record = next(parse(file, 'fasta'))

sequence_record_dict = to_dict(parse(file, 'fasta'))

#convert FIRST (0) dictionary keys to a list 
sequence_record = sequence_record_dict[list(sequence_record_dict)[0]]

#print(sequence_record.id)

rna_sequence = transcribe(sequence_record.seq)

protein_sequence = rna_sequence.translate()

#print(protein_sequence)

blastn_store = NCBIWWW.qblast("blastn", "nt", sequence_record.seq)
blastn_results = blastn_store.read() 
print(blastn_results)





