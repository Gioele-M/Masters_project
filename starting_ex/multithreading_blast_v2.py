from Bio import SeqIO
from Bio.Blast import NCBIWWW
from itertools import repeat
import time
import concurrent.futures


def open_fasta(filename):
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    return sequence_record


def run_blast(fasta_record, taxid, query_size = 20):
    print(taxid)
    entrez = f'txid{taxid}[ORGN]'
    result_handler = NCBIWWW.qblast('blastn', 'nt', fasta_record.seq, entrez_query= entrez, hitlist_size=query_size)
    result_storer = result_handler.read()
    print(f'{taxid} done')
    return result_storer


if __name__ ==  '__main__':
    
    start = time.perf_counter()
    fasta_record = open_fasta('starting_ex/human_mx1.fas')
    #old_taxid_list = ['9592', '9527', '756884', '314147', '91950', '9544', '37011', '147650', '60711', '54600']

    taxid_list = ['9606', '9597', '9593', '9600', '9601', '61853', '9546', '9544', '9541', '54180']


    iterable = zip(repeat(fasta_record), taxid_list)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = [executor.submit(run_blast, fasta_record, taxid) for taxid in taxid_list]
    for f in concurrent.futures.as_completed(results):
        print(f'done!!')

    end = time.perf_counter()
    print(f'Finished in {round(end-start,2)}')



#perc_ident=