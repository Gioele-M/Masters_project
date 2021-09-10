
import multiprocessing
import os
import pandas as pd
import numpy as np
from Bio import SeqIO, SeqRecord, Seq, SearchIO, AlignIO, Phylo
from Bio.Blast import NCBIWWW, NCBIXML
import Bio.Entrez


def open_fasta(filename):
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    return sequence_record

def run_blast(fasta_record, taxid):
    print(taxid)
    entrez = f'txid{taxid}[ORGN]'
    result_handler = NCBIWWW.qblast('blastn', 'nt', fasta_record, entrez_query= entrez)
    result_storer = result_handler.read()
    print(f'{taxid} done')
    return result_storer


if __name__ ==  '__main__':
    
    fasta_record = open_fasta('starting_ex/human_mx1.fas')

    taxid_list = ['9592','9527', '40674', '314147']


    threads=2
    pool = multiprocessing.Pool(threads)
    results = pool.starmap(run_blast,[(fasta_record.seq, taxid_list[0]), (fasta_record.seq, taxid_list[1])]) #!!!!!!!
    pool.close()
    pool.join()









#file = open('human_mx1.fas') 

#get the first itterable with next() 
#sequence_record = next(parse(file, 'fasta'))

#sequence_record_dict = to_dict(parse(file, 'fasta'))

#convert FIRST (0) dictionary keys to a list 
#sequence_record = sequence_record_dict[list(sequence_record_dict)[0]]

#print(sequence_record.id)

#rna_sequence = transcribe(sequence_record.seq)

#protein_sequence = rna_sequence.translate()

#print(protein_sequence)

'''blastn_store = NCBIWWW.qblast("blastn", "nt", fasta_record.seq)
blastn_results = blastn_store.read() 
print(blastn_results)'''


