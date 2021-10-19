import multiprocessing
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from itertools import repeat
import time
import datetime



def open_fasta(filename):
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    return sequence_record


def run_blast(fasta_record, taxid, query_size = 20):
    #print(taxid)
    entrez = f'txid{taxid}[ORGN]'
    result_handler = NCBIWWW.qblast('blastn', 'nt', fasta_record.seq, entrez_query= entrez, hitlist_size=query_size)
    result_storer = result_handler.read()
    #print(f'{taxid} done')
    return result_storer


if __name__ ==  '__main__':
    
    start = time.perf_counter()

    fasta_record = open_fasta('human_mx1.fas')
    taxid_list = ['9606', '9597', '9593', '9600', '9601', '61853', '9546', '9544', '9541', '54180']

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

    with open('test_multiprocessing.xml', 'w') as handle:
        handle.write(storer)

    end = time.perf_counter()
    timestamp = datetime.datetime.now()

    print(f'Multiprocessing of {len(taxid_list)} taxids finished in {round(end-start,2)}; Time {timestamp}')