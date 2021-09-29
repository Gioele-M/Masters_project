import os
from sys import stderr, stdout
import tempfile
import subprocess
from warnings import catch_warnings
import pandas as pd
import numpy as np
from Bio import SeqIO, SeqRecord, Seq, SearchIO, AlignIO, Phylo
from Bio.Blast import NCBIWWW, NCBIXML
import Bio.Entrez
from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from pandas.core import construction
from pandas.core.frame import DataFrame

import re
import time


#Function to open fasta file of imput
def open_fasta(filename) -> SeqRecord:
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    return sequence_record


#Function to run BLAST with taxid list
def blastp_with_list(sequence, list_taxid = [], query_size = 200):
    result_handler, result_storer = None, None
    #If list is empty run query without specific taxid
    if len(list_taxid) <1:
        result_handler = NCBIWWW.qblast('blastp', 'nr', sequence, hitlist_size=query_size) 
        result_storer = result_handler.read()
    #Prepare string of Entrez and parse it to qblast
    else:
        entrez_query = ''
        for taxid in list_taxid:
            entrez_query += f'txid{taxid}[ORGN]'
            if taxid != list_taxid[-1]:
                entrez_query += ' OR '
        result_handler = NCBIWWW.qblast('blastp', 'nr', sequence, entrez_query= entrez_query, hitlist_size=query_size)
        result_storer = result_handler.read()
    return result_storer



#Single blast query for threading
#Threading to be implemented 
def blastp_single_taxid(sequence, taxid, query_size = 20):
    entrez_query = f'txid{taxid}[ORGN]'
    result_handler = NCBIWWW.qblast('blastp', 'nr', sequence, entrez_query= entrez_query, hitlist_size=query_size)    
    result_storer = result_handler.read()
    return result_storer


#Function to create a temporary file from string and read it with SearchIO.read
def xml_string_to_handler(string):
    #make temporary file
    tmp = tempfile.NamedTemporaryFile(mode='a+')
    #write string
    tmp.write(string)
    handler = SearchIO.read(tmp.name, 'blast-xml')
    tmp.close()
    return handler


#Creation of a dictionary with all HSPS
#Bitscore, Bitscore raw and Evalue appear twice given that are difficult to average. For later QC only the values of the first hsp are considered
def blast_to_dictionary(blastresult):
    blast_dictionary = {'ID' : [], 'Description' : [], 'Seq_length' : [], 'Accession' : [], 'Bitscore' : [], 'Evalue' : [], 'Tot_aln_span':[], 'Identity' :[]}
    #Loop through results 
    for result in blastresult:
        blast_dictionary['ID'].append(result.id)
        blast_dictionary['Description'].append(result.description)
        blast_dictionary['Seq_length'].append(result.seq_len)
        blast_dictionary['Accession'].append(result.accession)
        #Store results of first HSP
        first_hsp = result.hsps[0]
        blast_dictionary['Bitscore'].append(first_hsp.bitscore)
        blast_dictionary['Evalue'].append(first_hsp.evalue)
        #Create variables to store results of multiple hsps
        all_alnspan, all_gapnum = [],[] 
        #Collect data of all hsps for each hit #REMOVE!!!!!!!!!
        for hsp in result.hsps:
            all_alnspan.append(int(hsp.aln_span))
            all_gapnum.append(int(hsp.gap_num))
        #Calculate total alignment span and gaps to calculate identity
        tot_alnspan, tot_gapnum = int(), int()
        seq_len = int(result.seq_len) #DOUBLE CHECK 
        for span in all_alnspan:
            tot_alnspan += span
        for gap in all_gapnum:
            tot_gapnum += gap
        identity = (tot_alnspan - tot_gapnum)/seq_len*100
        blast_dictionary['Tot_aln_span'].append(tot_alnspan)
        blast_dictionary['Identity'].append(round(identity, 3))
    return blast_dictionary


#Filter DF based on Evalue, sequence length and identity
def filter_df_blast(df:DataFrame, query:SeqRecord, evalue = 10**-10, difference_from_query = 50, identity_threshold=50) -> DataFrame:
    remove_index = []
    #Check if main HSP is significant for threshold, store indexes of non-significant hits
    for i in range(len(df.index)):
        eval = df.iloc[i]['Evalue']
        if float(eval) > evalue:
            remove_index.append(str(i))
        
    #Copy DF
    df_to_return = df    
    #Remove indexes if list is not empty
    if len(remove_index)>0:
        df_to_return = df_to_return.drop(df.index[int(remove_index)])
        df_to_return = df_to_return.reset_index(drop=True)
    
    #Filter by +/- % of query length
    #Obtain query length
    query_length = len(query.seq)
    #Determine upper/lower threshold of acceptance for sequences 
    lower, upper = query_length*(difference_from_query/100), query_length*((difference_from_query + 100)/100)
    #Filter DF
    df_to_return = df_to_return[df_to_return['Seq_length'] > lower]
    df_to_return = df_to_return.reset_index(drop=True)
    df_to_return = df_to_return[df_to_return['Seq_length'] < upper]
    df_to_return = df_to_return.reset_index(drop=True)

    #Filter by identity
    df_to_return = df_to_return[df_to_return['Identity'] > identity_threshold]
    df_to_return = df_to_return.reset_index(drop=True)
    
    return df_to_return


#Retrieve all result's entries from efetch
def retrieve_all_efetch(list_of_entries, email):
    Bio.Entrez.email = email
    list_of_sequences = []
    for entry in list_of_entries:
        handler = Bio.Entrez.efetch(db='protein', id=entry, rettype = 'fasta',retmode = 'xml', retmax=1) #Returns JSON regardless
        gb_info = Bio.Entrez.read(handler, 'text')#Returns nested lists and dictionaries 
        list_of_sequences.append(gb_info)
    return list_of_sequences



#Efetch to dictionary parser
def efetch_protein_to_dictionary(list_of_efetch):
    #Declare new dictionary
    dictionary = {'Accession':[],'Protein_ID':[], 'Taxid':[], 'Organism_name':[], 'Description':[], 'Seq_length':[], 'Prot_sequence':[]}
    for wrapper in list_of_efetch:
        #Cast into dictionary to avoid random exception
        result = dict(wrapper[0])
        acc_ver = result['TSeq_accver']
        accession = acc_ver.split('.')
        dictionary['Accession'].append(accession[0])
        dictionary['Protein_ID'].append(result['TSeq_accver'])
        dictionary['Taxid'].append(result['TSeq_taxid'])
        dictionary['Organism_name'].append(result['TSeq_orgname'])
        dictionary['Description'].append(result['TSeq_defline'])
        dictionary['Seq_length'].append(result['TSeq_length'])
        dictionary['Prot_sequence'].append(result['TSeq_sequence'])
        
    return dictionary


def filter_df_taxon(df:DataFrame, n_of_sequences = 1) -> DataFrame:
    #Make a list of all retrieved taxons
    retrieved_taxids = df['Taxid'].tolist() 
    #Make list from dictionary to eliminate duplicates 
    retrieved_taxids = list(dict.fromkeys(retrieved_taxids))
        
    #Declare empty DF
    df_toreturn = pd.DataFrame()

    #For each taxon create a DF, sort and get best result to append to df_toreturn
    for taxid in retrieved_taxids:
        #Create temporary DF exclusive to taxon
        temp_df = df[df['Taxid'] == taxid]
        temp_df = temp_df.sort_values(['Evalue', 'Identity', 'Bitscore'], ascending=[True, False, False]) #('Identity', ascending=False)
        if len(temp_df) > 0:
            if len(temp_df) > n_of_sequences:
                df_toreturn = df_toreturn.append(temp_df[:n_of_sequences])
            else:
                df_toreturn = df_toreturn.append(temp_df)

    #Refactor indexes
    df_toreturn = df_toreturn.reset_index(drop=True)

    return df_toreturn


#Prepare string for alingment file
def fasta_for_alignment(query:SeqRecord, df:DataFrame) -> str:
    #Initialise string and add 
    string = ''
    string += f">Query:{query.id}\n"
    string += f"{query.seq}\n"
    #Add all sequence IDs and aa sequence
    for i in range(len(df['Accession'])):
        #Check if protein is missing ------------------------------------------------Add to report? ERROR.LOG FILE!!!!!
        if len(df['Prot_sequence'][i]) <1:
            continue
        string += f">{df['ID'][i]}\n" # {df['Organism_name'][i]} {df['Description'][i]}
        string += f"{df['Prot_sequence'][i]}"
        #Check if is not the last element in the list 
        if i != len(df['Accession']):
            string += '\n'
    
    return string


#Run mafft by saving the fasta sequences for alignment in a file and passing it to mafft 
def run_mafft_saving_file(fasta:str, mafft_directory:str, filename:str) -> str:
    #Write fasta file for alignment
    file = f'{filename}.fasta'
    with open(file, 'w') as handle:
        handle.write(fasta)
    #Parse and run with mafft 
    command_list = [mafft_directory, '--distout', f'{file}']
    process = subprocess.Popen(command_list, universal_newlines= True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()

    #Save standard error of mafft
    with open('aligned_mafft.stderr', 'w') as handle:
        handle.write(stderr)
    #Retrun alignment 
    return stdout


#Function to get protein sequence and save fasta file given an accession number
#The function recycles functions used with lists 
def get_fasta_from_accession(accession, output_name):
    #Cast into list to recycle code 
    accession_query_list = [accession]
    #Get the efetch handler 
    query_protein_efetch = retrieve_all_efetch(accession_query_list, 'A.N.Other@example.com')
    #Make a dictionary
    dictionary_query = efetch_protein_to_dictionary(query_protein_efetch)

    #Make and save fasta file 
    fasta_string_query = f">{dictionary_query['Accession'][0]} \n {dictionary_query['Prot_sequence'][0]}"

    fasta_file_name = output_name + '_sequence.fasta'
    with open(fasta_file_name, 'w') as handle:
        handle.write(fasta_string_query)

    fasta_record_q = open_fasta(fasta_file_name)
    return fasta_record_q


def tree_from_alignment(alignment):
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    return tree





#MAIN 
if __name__ == '__main__':
    #USER INPUTS
    filename = 'sry_protein.fasta' 
    taxid_list = [] #['9592']
    mafft_directory = r'/Users/Gioele/miniconda3/bin/mafft'
    email = 'A.N.Other@example.com'
    output_name = 'draft_v4_1'

    #OPTIONAL USER INPUTS
    #type_of_query_list = ['nucleotide', 'amino_acid', 'accession'] #Default amino_acid
    type_of_query = 'amino_acid'                        #DETERMINE YOURSELF!!!!!! CHECK IF NOT A FASTA FILE, THEN CHECK FOR AA IN SEQ
                                                        #CAN CHECK IF THERE IS A .FAS/.FASTA IN THE NAME +++ ERROR LOG FILE 
    #Used if accession option is on 
    accession = 'QBA69874'


    local_query = False
    threading = False
    query_size = 1000

    evalue_threshold = 10**-10
    len_threshold = 50
    identity_threshold = 50

    sequences_per_taxon = 1


    #Filenames:
    #Implemented
    output_df = f'{output_name}_df.csv'
    output_alignment = f'{output_name}_alignment.fasta'
    output_xml_tree = f'{output_name}_xml_tree.xml'
    #To add
    output_tree_newick = f'{output_name}_newick_tree.nwk' #!!!!
    output_tree_jpg = f'{output_name}_tree_image.jpg'








    #START OF THE SCRIPT
    start = time.perf_counter()

    #Declare fasta record, check loop to fill it 
    fasta_record = None
    #Open fasta record appropriately 
    if type_of_query == 'nucleotide':
        #Get protein from nucleotide 
        fasta_nucleotide = open_fasta(filename) #filename 
        #Further QC? 
        fasta_record = fasta_nucleotide.translate(to_stop=True)
    elif type_of_query == 'accession':
        #Get fasta record from accession 
        fasta_record = get_fasta_from_accession(accession, output_name)
    else: 
        fasta_record = open_fasta(filename)

    
    print('Opened record')
    print(f'In {round(time.perf_counter()-start,2)}')
    

    #NEED TO DOUBLECHECK AND IMPLEMENT BLAST WITH THREADING 

    #Blast with list 
    blast_results = blastp_with_list(fasta_record.seq, query_size=query_size, list_taxid=taxid_list)
    
    #Parse handler to xml_string_to_handler to be able to feed it in SearchIO.read
    handler_blast = xml_string_to_handler(blast_results)

    #Make dictionary from handler
    dictionary_blast = blast_to_dictionary(handler_blast)

    #Transfor dictionary into DF
    results_df_blast = pd.DataFrame.from_dict(dictionary_blast)

    print('Prepared DF BLAST')
    print(f'In {round(time.perf_counter()-start,2)}')


    #Filter DF by E-value / Bitscore / Identity
    filtered_blast_df = filter_df_blast(results_df_blast, fasta_record, evalue=evalue_threshold, difference_from_query=len_threshold, identity_threshold=identity_threshold)
 
    #Retrieve all results with efetch 
    #Retrieve accession list
    accession_list = filtered_blast_df['Accession'].tolist()
    #Get results 
    efetch_result = retrieve_all_efetch(accession_list, email)

    # Make dictionary from efetch fasta results
    dictionary_efetch = efetch_protein_to_dictionary(efetch_result)

    #Transform dictionary into DF
    efetch_df = pd.DataFrame.from_dict(dictionary_efetch)

    print('Done Efetch')
    print(f'In {round(time.perf_counter()-start,2)}')



    #Merge blast and efetch DataFrames
    left = filtered_blast_df.loc[:,['Accession', 'ID','Seq_length', 'Evalue', 'Bitscore', 'Tot_aln_span', 'Identity']]
    right = efetch_df.loc[:,['Accession', 'Taxid', 'Organism_name', 'Description', 'Prot_sequence']]
    combined_df = pd.merge(left, right, on='Accession')


    #Filter DF based on taxons
    filtered_df = filter_df_taxon(combined_df, n_of_sequences=sequences_per_taxon)

    #Save DF of sequences that are going to be aligned 
    '''with open(output_df, 'w') as handle:
        handle.write(filtered_df.to_csv())'''
    filtered_df.to_csv(output_df, index = False)
    
    #####df = df[ ['Mid'] + [ col for col in df.columns if col != 'Mid' ] ]
    ###############
    tsv_df = output_df[['ID'] + [col for col in output_df.columns if col!= 'ID']]
    #Add TSV!!! Rearrange columns for ProteinID first!!!!!
    filtered_df.to_csv(tsv_df, sep = '\t', index = False)


    #Prepare for multiple alingment
    fasta_string_for_alignment = fasta_for_alignment(fasta_record, filtered_df)
    
    #Run mafft alignment saving file
    mafft_alignment = run_mafft_saving_file(fasta_string_for_alignment, mafft_directory, 'multiple_seq_fasta')

    #Save mafft alignment
    with open(output_alignment, 'w') as savefile:
        savefile.write(mafft_alignment)

    print('Alignment done')
    print(f'In {round(time.perf_counter()-start,2)}')


    #Open AlignIO from fasta file 
    alignment = AlignIO.read(output_alignment, 'fasta')

    #Construct tree from alignment 
    #To IMPLEMENT DIFFERENT METHODS FOR TREE BUILDING 
    tree = tree_from_alignment(alignment)  


    #Write tree in xml
    Phylo.write(tree, output_xml_tree, 'phyloxml')

    Phylo.write(tree, output_tree_newick, 'newick')


    print('Done!')  #MORE DETAILED!!!!! SO MANY SEQUENCES FOUND FOR THIS MANY TAXA, PRINTED IN THIS FILE, THE LOG IS THAT FILE ETC.....

    end = time.perf_counter()
    print(f'Finished in {round(end-start,2)}')












    '''
            
        #make sure it's a coding sequence
        if len(seq)%3 != 0:
            print("\n\n\\(*'o'*)/\tOops! Sequence %s has length not a multiple of 3...\n\n"%str(rec.id))
            return
        
        
        #save amino acid sequence as string
        aa = str(Seq(seq).translate())
        
        #remove stop codons at the end
        if aa[-1] == '*':
            aa = aa[:-1]
            seq = seq[:-3]
        
        #check for internal stop codons
        if '*' in aa:
            print("\n\n\\(*'o'*)/\tOops! Sequence %s has internal stop codons...\n\n"%str(rec.id))
            return  


    
    '''
