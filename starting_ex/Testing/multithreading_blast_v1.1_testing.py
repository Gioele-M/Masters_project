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
    taxid_list = [] #['9606', '9597', '9593', '9600', '9601', '61853', '9546', '9544', '9541', '54180']
    new_list = [163843, 109228, 98314, 98318, 109229, 163861, 196630, 1081385, 1704005, 65606, 65607, 65609, 65610, 65613, 65614, 65615, 65616, 65618, 65621, 65624, 65626, 65628, 65629, 65631, 65632, 65634, 1851651, 458857, 458858, 458859, 458860, 458861, 458862, 1048687, 1048688, 655476, 1649342, 2750743, 1851652, 2588810, 92866, 1392492, 65694, 1851653, 190511, 2162857, 2162858, 2162859, 2162860, 2162861, 2162862, 2162864, 2162865, 2162866, 2162901, 2162903, 2162904, 1147103, 1147105, 1147106, 655592, 491753, 557290, 491755, 1048812, 1048813, 1294088, 48018, 327948, 327957, 1966371, 1966372, 2588968, 2588969, 1605932, 1966381, 1966382, 1966383, 1966384, 1966385, 1513475, 131388, 262468, 491852, 491853, 1048917, 1048918, 1048919, 1048920, 1048921, 1048922, 1048923, 1081694, 491872, 491874, 491876, 491878, 1343853, 1343854, 491888]
    for x in new_list:
        taxid_list.append(str(x))

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
        try:
            tmp = tempfile.NamedTemporaryFile(mode='a+') 
            tmp.write((result))
            results_handler = SearchIO.read(tmp.name, 'blast-xml')
            list_of_handlers.append(results_handler)
            tmp.close()
        except Exception:
            print('A result couldnt be parsed')


    storer = ''
    
    for request in results:
        storer += str(request)
        storer += '\n'

    with open('test_threading.xml', 'w') as handle:
        handle.write(storer)

    end = time.perf_counter()
    print(f'Threading of {len(taxid_list)} taxids finished in {round(end-start,2)}')

