import threading
from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW
from itertools import repeat
import time
import tempfile
import pandas as pd


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

#Creation of dictionary with all HSPs adding identity and other metrics

#Hits without at least a significant HSP are excluded

def blast_to_dictionary_plus_metrics(blastresult):
    blast_dictionary = {'ID' : [], 'Description' : [], 'Seq_length' : [], 'Accession' : [], 'Bitscore' : [], 'Bitscore_raw' : [], 
    'Evalue' : [], 'Hit_start' : [], 'Hit_end' : [], 'Query_frame' : [], 'Gap_num' : [], 'Aln_span' : [], 'Tot_aln_span':[], 'Identity' :[]}
    for result in blastresult:
        blast_dictionary['ID'].append(result.id)
        blast_dictionary['Description'].append(result.description)
        blast_dictionary['Seq_length'].append(result.seq_len)
        blast_dictionary['Accession'].append(result.accession)
        bitscore, bitscore_raw, evalue, hitstart, hitend, queryframe, gapnum, alnspan = '','','','','','','',''
        all_alnspan, all_gapnum = [],[] #seq_len is not required
        for hsp in result.hsps:
            bitscore += str(hsp.bitscore)
            bitscore_raw += str(hsp.bitscore_raw)
            evalue += str(hsp.evalue)
            hitstart += str(hsp.hit_start)
            hitend += str(hsp.hit_end)
            queryframe += str(hsp.query_frame)
            gapnum += str(hsp.gap_num)
            alnspan += str(hsp.aln_span)
            if hsp != result.hsps[-1]:
                bitscore += '/'
                bitscore_raw += '/' 
                evalue += '/'
                hitstart += '/'
                hitend += '/'
                queryframe += '/'
                gapnum += '/' 
                alnspan += '/' #I know it's not neat but it would give an error otherwise :(
            all_alnspan.append(int(hsp.aln_span))
            all_gapnum.append(int(hsp.gap_num))
        blast_dictionary['Bitscore'].append(bitscore)
        blast_dictionary['Bitscore_raw'].append(bitscore_raw)
        blast_dictionary['Evalue'].append(evalue)
        blast_dictionary['Hit_start'].append(hitstart)
        blast_dictionary['Hit_end'].append(hitend)
        blast_dictionary['Query_frame'].append(queryframe)
        blast_dictionary['Gap_num'].append(gapnum)
        blast_dictionary['Aln_span'].append(alnspan)
        tot_alnspan, tot_gapnum = int(), int()
        seq_len = int(result.seq_len)
        for span in all_alnspan:
            tot_alnspan += span
        for gap in all_gapnum:
            tot_gapnum += gap
        identity = (tot_alnspan - tot_gapnum)/seq_len*100
        blast_dictionary['Tot_aln_span'].append(tot_alnspan)
        blast_dictionary['Identity'].append(round(identity, 3))
    return blast_dictionary


if __name__ ==  '__main__':
    
    start = time.perf_counter()
    fasta_record = open_fasta('starting_ex/human_mx1.fas')
    #old_taxid_list = ['9592', '9527', '756884', '314147', '91950', '9544', '37011', '147650', '60711', '54600']
    taxid_list = ['9606', '9597', '9593', '9600', '9601', '61853', '9546', '9544', '9541', '54180']

    #iterable = zip(repeat(fasta_record), taxid_list)

    threads = []
    results = []
    list_of_handlers = []
    data_frame = pd.DataFrame()

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


    for handler in list_of_handlers:
        dictionary = blast_to_dictionary_plus_metrics(handler)
        data_frame.append(pd.DataFrame.from_dict(dictionary))
    








    storer = ''
    
    for request in results:
        storer += str(request)
        storer += '\n'

    with open('test_threading.xml', 'w') as handle:
        handle.write(storer)

    end = time.perf_counter()
    print(f'Finished in {round(end-start,2)}')






#arg_line = " ".join(sys.argv[1:])
'''
def parse(arg_line: str) -> Dict[str, str]:
    args: Dict[str, str] = {}
    if match_object := args_pattern.match(arg_line):
        args = {k: v for k, v in match_object.groupdict().items()
                if v is not None}
    return args
'''
#TO TRY!!!!!!