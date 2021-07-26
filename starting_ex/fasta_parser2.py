import pandas as pd
import numpy as np
from Bio import SeqIO, SeqRecord, Seq, SearchIO
from Bio.Blast import NCBIWWW

def mk_sequence_record(file, format='fasta'):
    seq_record_dictionary = SeqIO.to_dict(SeqIO.parse(file, format)) #maybe best with read instead of parse
    sequence_record = seq_record_dictionary[list(seq_record_dictionary)[0]]
    return sequence_record


def covert_to_protein(seq):
    rna = Seq.transcribe(seq.seq)
    protein = rna.translate(to_stop=True)
    return protein #returns Seq object


def blast_search_xml_save(sequence, filename, taxid = '',type='blastn', database = 'nt'):
    if len(taxid) >0:
        handler = NCBIWWW.qblast(type, database, sequence.seq, entrez_query=taxid)
    else:
        handler = NCBIWWW.qblast(type, database, sequence.seq)
    storer = handler.read()
    with open(f'{filename}.xml', 'w') as savefile:
        savefile.write(storer)



def single_result_xml_read_to_pd_dict(filename):
    result = SearchIO.read(f'{filename}.xml', 'blast-xml')
    dictionary = SearchIO.to_dict(result)
    df = pd.DataFrame.from_dict(dictionary)
    return df



seq = mk_sequence_record("human_mx1.fas")
blast_search_xml_save(seq, 'try1')
dataframe = single_result_xml_read_to_pd_dict('try1')
print(dataframe)