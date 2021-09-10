import multiprocessing
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from itertools import repeat
import time



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
    taxid_list = ['9592', '9527', '40674', '314147', '9531', '9544', '2008792']

    threads=multiprocessing.cpu_count()-1
    pool = multiprocessing.Pool(threads)
    iterable = zip(repeat(fasta_record), taxid_list)
    results = pool.starmap(run_blast, iterable) #!!!!
    pool.close()
    pool.join()

    storer = ''
    
    for request in results:
        storer += str(request)
        storer += '\n'

    with open('test_multithread.xml', 'w') as handle:
        handle.write(storer)

    end = time.perf_counter()
    print(f'Finished in {round(end-start,2)}')