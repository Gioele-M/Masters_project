import multiprocessing
from Bio import SeqIO
from Bio.Blast import NCBIWWW




def open_fasta(filename):
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    return sequence_record


def blast_single_taxid_parallel(sequence, filename, query_size = 50, list_taxid = []):
    
    
    def run_blast(taxid):
        entrez = f'txid{taxid}[ORGN]'
        result_handler = NCBIWWW.qblast('blastn', 'nt', sequence, entrez_query= entrez, hitlist_size=query_size)
        result_storer = result_handler.read()
        print(f'{taxid} done')
        return result_storer
    
    
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count()-1)
    storer_list = pool_obj.map(run_blast, list_taxid) #for x in list_taxid
    return storer_list



if __name__ ==  '__main__':
    
    fasta_record = open_fasta('starting_ex/human_mx1.fas')

    taxid_list = ['9592','9527']#, '40674', '314147']

    results = blast_single_taxid_parallel(fasta_record.seq, 'useless_now', query_size=25, list_taxid=taxid_list)



    '''threads=2
    pool = multiprocessing.Pool(threads)
    results = pool.starmap(run_blast,[(fasta_record.seq, taxid_list[0]), (fasta_record.seq, taxid_list[1])]) #!!!!
    pool.close()
    pool.join()'''

    storer = ''
    
    for request in results:
        storer += str(request)
        storer += '\n'

    with open('test_multithread_v3.xml', 'w') as handle:
        handle.write(storer)