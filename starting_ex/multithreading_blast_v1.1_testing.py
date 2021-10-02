import threading
from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW
from itertools import repeat
import time
import tempfile

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
    fasta_record = open_fasta('human_mx1.fas')
    taxid_list = ['9606', '9597', '9593', '9600', '9601', '61853', '9546', '9544', '9541', '54180']

    #iterable = zip(repeat(fasta_record), taxid_list)

    threads = []
    results = []
    list_of_handlers = []

    for taxid in taxid_list:
        t = threading.Thread(target=lambda: results.append(run_blast(fasta_record, taxid))) #Number of threads
        t.start()
        threads.append(t)

    for thread in threads:
        thread.join()



    #Make temporary file, write string, pass file to SearchIO.read, pass the handler to list, close temporary file
    for result in results:
        tmp = tempfile.NamedTemporaryFile(mode='a+') 
        tmp.write((result))
        results_handler = SearchIO.read(tmp.name, 'blast-xml')
        list_of_handlers.append(results_handler)
        tmp.close()


    storer = ''
    
    for request in results:
        storer += str(request)
        storer += '\n'

    with open('test_threading.xml', 'w') as handle:
        handle.write(storer)

    end = time.perf_counter()
    print(f'Threading of {len(taxid_list)} taxids finished in {round(end-start,2)}')

