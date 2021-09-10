import threading
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
    #old_taxid_list = ['9592', '9527', '756884', '314147', '91950', '9544', '37011', '147650', '60711', '54600']
    taxid_list = ['9606']#, '9597', '9593', '9600', '9601', '61853', '9546', '9544', '9541', '54180']

    iterable = zip(repeat(fasta_record), taxid_list)

    threads = []

    '''for i in range(len(taxid_list)):
        t = threading.Thread(target=run_blast, args=[fasta_record, taxid_list[i]])
        t.start()
        threads.append(t)'''
    #WORKS TOO!!!!


    for i in iterable:
        t = threading.Thread(target=run_blast, args=i)
        t.start()
        threads.append(t)

    for thread in threads:
        thread.join()


    '''for thread in threads:
        print(thread)'''

    storer = ''
    
    for request in thread:
        storer += str(request)
        storer += '\n'

    with open('test_threading.xml', 'w') as handle:
        handle.write(storer)

    end = time.perf_counter()
    print(f'Finished in {round(end-start,2)}')