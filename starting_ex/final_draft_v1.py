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

import time


 
#Function to open fasta file of imput
def open_fasta(filename):
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    return sequence_record


#Function to run BLAST with taxid list
def blastn_with_list(sequence, list_taxid = [], query_size = 200, blast_type = 'blastn'):
    result_handler, result_storer = None, None
    #If list is empty run query without specific taxid
    if len(list_taxid) <1:
        result_handler = NCBIWWW.qblast('blastn', 'nt', sequence, hitlist_size=query_size) 
        result_storer = result_handler.read()
    #Prepare string of Entrez and parse it to qblast
    else:
        entrez_query = ''
        for taxid in list_taxid:
            entrez_query += f'txid{taxid}[ORGN]'
            if taxid != list_taxid[-1]:
                entrez_query += ' OR '
        result_handler = NCBIWWW.qblast('blastn', 'nt', sequence, entrez_query= entrez_query, hitlist_size=query_size)
        result_storer = result_handler.read()
    return result_storer


#Single blast query for threading
def blast_single_taxid(sequence, taxid, query_size = 20):
    entrez_query = f'txid{taxid}[ORGN]'
    result_handler = NCBIWWW.qblast('blastn', 'nt', sequence, entrez_query= entrez_query, hitlist_size=query_size)    
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
    blast_dictionary = {'ID' : [], 'Description' : [], 'Seq_length' : [], 'Accession' : [], 'Bitscore' : [], 'Bitscore_raw' : [], 
    'Evalue' : [], 'Bitscore_combined':[], 'Bitscore_raw_combined' : [], 'Evalue_combined':[],'Hit_start' : [], 'Hit_end' : [], 'Query_frame' : [], 'Gap_num' : [], 'Aln_span' : [], 'Tot_aln_span':[], 'Identity' :[]}
    #Loop through results 
    for result in blastresult:
        blast_dictionary['ID'].append(result.id)
        blast_dictionary['Description'].append(result.description)
        blast_dictionary['Seq_length'].append(result.seq_len)
        blast_dictionary['Accession'].append(result.accession)
        #Store results of first HSP
        first_hsp = result.hsps[0]
        blast_dictionary['Bitscore'].append(first_hsp.bitscore)
        blast_dictionary['Bitscore_raw'].append(first_hsp.bitscore_raw)
        blast_dictionary['Evalue'].append(first_hsp.evalue)
        #Create variables to store results of multiple hsps
        bitscore, bitscore_raw, evalue, hitstart, hitend, queryframe, gapnum, alnspan = '','','','','','','',''
        all_alnspan, all_gapnum = [],[] 
        #Collect data of all hsps for each hit
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
        blast_dictionary['Bitscore_combined'].append(bitscore)
        blast_dictionary['Bitscore_raw_combined'].append(bitscore_raw)
        blast_dictionary['Evalue_combined'].append(evalue)
        blast_dictionary['Hit_start'].append(hitstart)
        blast_dictionary['Hit_end'].append(hitend)
        blast_dictionary['Query_frame'].append(queryframe)
        blast_dictionary['Gap_num'].append(gapnum)
        blast_dictionary['Aln_span'].append(alnspan)
        #Calculate total alignment span and gaps to calculate identity
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
        

#Retrieve all result's entries from Genbank
def retrieve_all_genbank(list_of_entries, email):
    Bio.Entrez.email = email
    list_of_sequences = []
    for entry in list_of_entries:
        handler = Bio.Entrez.efetch(db='nucleotide', id=entry, rettype = 'gb', retmode = 'xml', retmax=1) #!!!Genbank format !! Returns JSON regardless
        #retmode = 'text'
        gb_info = Bio.Entrez.read(handler, 'genbank')
        list_of_sequences.append(gb_info)
    return list_of_sequences



#Create a dictionary from the genbank results
def genbank_to_dictionary(list_of_genbank):
    #Define function to retrieve values from the object. If the value doesn't exist an empty string is returned to avoid errors
    def find_value(result, parameter):
        to_return = ''
        try:
            to_return = result[parameter]
        except KeyError:
            print(f'Could not find {parameter}')
        return to_return


    
    #Define functions to find taxon, protein ID and protein Sequence in the JSON nest
    #Qualifiers names for the attributes are not flexible
    def find_taxon(list_feature_table):
        taxon = ''
        #Loop through all the features in the list of features
        for feature in list_feature_table:
            #If the feature is source hold the wrapper and loop through it
            if feature['GBFeature_key'] == 'source':
                source_wrap = feature['GBFeature_quals']
                for qualifier in source_wrap:
                    #Check if the qualifier is actually a taxon and store it in the variable
                    if qualifier['GBQualifier_name'] == 'db_xref' and qualifier['GBQualifier_value'].startswith('taxon'):
                        taxon += (qualifier['GBQualifier_value'])
        #Comment: artificial constructs have often more than 2 taxons, this will be removed later
        return taxon    

    def find_prot_id(list_feat_table):
        prot_id = ''
        for feature in list_feat_table:
            if feature['GBFeature_key'] == 'CDS':
                wrapper = feature['GBFeature_quals']
                for qualifier in wrapper:
                    if qualifier['GBQualifier_name'] == 'protein_id':
                        prot_id += qualifier['GBQualifier_value']
        return prot_id

    def find_prot_seq(list_feat_table):
        prot_seq = ''
        for feature in list_feat_table:
            if feature['GBFeature_key'] == 'CDS':
                wrapper = feature['GBFeature_quals']
                for qualifier in wrapper:
                    if qualifier['GBQualifier_name'] == 'translation':
                        prot_seq += qualifier['GBQualifier_value']
        return prot_seq

    #Declare new dictionary
    dictionary_gen = {'Accession' : [], 'Accession_version': [], 'Gene_length' : [], 'Strandedness': [], 'Molecule_type':[], 'Organism':[], 'Taxonomy':[],
    'Nuc_sequence':[], 'Taxon':[], 'Protein_ID':[], 'Prot_sequence':[]} #'N_of_references':[]
    #Loop through genbank results and fill dictionary 
    
    #Pass all values through a function that returns 
    for gen in list_of_genbank:
        first_wrapper = gen[0]
        dictionary_gen['Accession'].append(find_value(first_wrapper, 'GBSeq_primary-accession'))
        dictionary_gen['Accession_version'].append(find_value(first_wrapper, 'GBSeq_accession-version'))
        dictionary_gen['Gene_length'].append(find_value(first_wrapper, 'GBSeq_length'))
        dictionary_gen['Strandedness'].append(find_value(first_wrapper, 'GBSeq_strandedness')) 
        dictionary_gen['Molecule_type'].append(find_value(first_wrapper, 'GBSeq_moltype'))
        dictionary_gen['Organism'].append(find_value(first_wrapper, 'GBSeq_organism'))
        dictionary_gen['Taxonomy'].append(find_value(first_wrapper, 'GBSeq_taxonomy'))
        dictionary_gen['Nuc_sequence'].append(find_value(first_wrapper, 'GBSeq_sequence'))

        #Old extraction raised error
        '''
        dictionary_gen['Accession'].append(first_wrapper['GBSeq_primary-accession'])
        dictionary_gen['Accession_version'].append(first_wrapper['GBSeq_accession-version'])
        dictionary_gen['Gene_length'].append(first_wrapper['GBSeq_length'])
        dictionary_gen['Strandedness'].append(first_wrapper['GBSeq_strandedness']) 
        dictionary_gen['Molecule_type'].append(first_wrapper['GBSeq_moltype'])
        dictionary_gen['Organism'].append(first_wrapper['GBSeq_organism'])
        dictionary_gen['Taxonomy'].append(first_wrapper['GBSeq_taxonomy'])
        dictionary_gen['Nuc_sequence'].append(first_wrapper['GBSeq_sequence'])
        '''

        dictionary_gen['Taxon'].append(find_taxon(first_wrapper['GBSeq_feature-table']))
        dictionary_gen['Protein_ID'].append(find_prot_id(first_wrapper['GBSeq_feature-table']))
        dictionary_gen['Prot_sequence'].append(find_prot_seq(first_wrapper['GBSeq_feature-table']))
    return dictionary_gen


def filter_df_taxon(df:DataFrame, n_of_sequences = 1) -> DataFrame:
    #Make a list of all retrieved taxons
    retrieved_taxons = df['Taxon'].tolist() 
    #Make list from dictionary to eliminate duplicates 
    retrieved_taxons = list(dict.fromkeys(retrieved_taxons))
    
    #Clean format
    retrieved_taxons = [taxid.replace('taxon:', '', 1) for taxid in retrieved_taxons] #Only remove first instance of taxon to be able to recognise synthetic constructs 

    #Declare final list to loop through
    final_list = []
    #Check if it has different taxids, if it does it's likely to be a synthetic construct so can be excluded
    for taxid in retrieved_taxons:
        #If taxid is 'clean' can be appended to the list
        if taxid.find('taxon') == -1:
            final_list.append(taxid)
        else:
            continue
    
    #Declare empty DF
    df_toreturn = pd.DataFrame()

    #For each taxon create a DF, sort and get best result to append to df_toreturn
    for taxon in final_list:
        #Create temporary DF exclusive to taxon
        temp_df = df[df['Taxon'] == taxon]
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
        #Check if protein is missing ------------------------------------------------Add to report?
        if len(df['Prot_sequence'][i]) <1:
            continue
        string += f">{df['Accession'][i]} {df['Organism'][i]}\n"
        string += f"{df['Prot_sequence'][i]}"
        #Check if is not the last element in the list 
        if i != len(df['Accession']):
            string += '\n'
    
    return string


def run_mafft_from_string(fasta:str, mafft_directory:str) -> str: 
    #Create temporary file to pass to mafft
    tmp = tempfile.NamedTemporaryFile(mode='a+')
    tmp.write(fasta)
    #Run mafft with subprocess
    command_list = [mafft_directory, '--distout', f'{tmp.name}']
    process = subprocess.Popen(command_list, universal_newlines= True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()

    #Save standard error of mafft
    with open('aligned_mafft_v1.stderr', 'w') as handle:
        handle.write(stderr)
    
    #Close temporary file
    tmp.close()
    #Retrun alignment 
    return stdout


def run_mafft_saving_file(fasta:str, mafft_directory:str, filename:str) -> str:
    file = f'{filename}.fasta'
    with open(file, 'w') as handle:
        handle.write(fasta)
    
    command_list = [mafft_directory, '--distout', f'{file}']
    process = subprocess.Popen(command_list, universal_newlines= True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()

    #Save standard error of mafft
    with open('aligned_mafft_v1.stderr', 'w') as handle:
        handle.write(stderr)

    #Retrun alignment 
    return stdout



#For some reason it throws an error depending on the sequences parsed

'''
#Create an alignment file from the alignment fasta string
def alignio_from_string(string:str):
    #Create temporary file to pass to AlignIO
    tmp = tempfile.NamedTemporaryFile(mode='a+', suffix='.fasta')
    tmp.write(string)
    #Open alignment
    filename = f'{tmp.name}'
    print(filename)
    alignment = AlignIO.read(open(filename), 'fasta')
    tmp.close()
    return alignment
'''



def tree_from_alignment(alignment):
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    return tree





if __name__ == '__main__':

    start = time.perf_counter()




    #User variables
    '''
    -Filename
    -Taxid List
    -EMAIL FOR GENBANK
    -MAFFT EXECUTABLE DIRECTORY 
    -Blast configuration
        -Local or cloud blast
        -Treading or not
            -!!!Number of threads !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Queue class
        -Query-size
        -Server used (blastn-nt/blastp-)
        -PARAMETERS FOR RESEARCHES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    -QC
        -Filter by Evalue
        -Filter by +-% query
        -Filter by identity
    -Further filtering
        -N of sequences per taxon

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    -If given mRNA is ambiguous/has different table, allow importing .fasta for AA seq

    '''
    
    #USER INPUTS
    filename = 'sry_gene.fasta' #'human_mx1.fas'
    taxid_list = [] #['9592']
    mafft_directory = r'/Users/Gioele/miniconda3/bin/mafft'
    email = 'A.N.Other@example.com'
    output_name = 'run1'
    

    #OPTIONAL USER INPUTS
    local_query = False
    threading = False
    query_size = 100
    server, database = '',''#FOR BLAST

    evalue_threshold = 10**-10
    len_threshold = 50
    identity_threshold = 50

    sequences_per_taxon = 3
    import_aa_sequence = True
    query_protein = 'sry_protein.fasta'
    
    






    #Implemented
    output_df = f'{output_name}_df.csv'
    output_alignment = f'{output_name}_alignment.fasta'
    output_xml_tree = f'{output_name}_xml_tree.xml'
    #To add
    output_tree_newick = f'{output_name}_newick_tree.newick'
    output_tree_jpg = f'{output_name}_tree_image.jpg'
    
    
    
    #Start of script


    handlers = []


    results_df = pd.DataFrame()

    #Open fasta record
    fasta_record = open_fasta(filename)
    
    print('Opened record')
    print(f'In {round(time.perf_counter()-start,2)}')


    #This has to be fixed 
    ''' 
    #Run blast, either in multi-list or single-threading - Append handlers to the list
    if not threading:
        #Taxid with list
        blast_results = blastn_with_list(fasta_record.seq, query_size=query_size, list_taxid=taxid_list)
        #Parse handler to xml_string_to_handler to be able to feed it in SearchIO.read
        handlers.append(xml_string_to_handler(blast_results))
    else:
        #Make list of threads and results
        threads, results = [], []
        #Create threads and join them
        for taxid in taxid_list:
            #Lambda function bypasses Thread functionality to return results from the function
            t = threading.Thread(target= lambda: results.append(blast_single_taxid(fasta_record.seq, taxid)))
            t.start()
            threads.append(t)
        for thread in threads:
            thread.join()
        #Parse handler to xml_string_to_handler to be able to feed it in SearchIO.read
        for result in results:
            handlers.append(xml_string_to_handler(result))


    print('Done Blast')

    #At this point all the results are in the list handlers
    #All the results will be converted to dictionary, then to DataFrame, and appended to the general DF
    #This solution was used to cover the cases for both inputs
    for handler in handlers:
        #Convert SearchIO handler to dictionary
        dictionary = blast_to_dictionary(handler)
        #Convert to DF and append it to general DF
        temp_df = pd.DataFrame.from_dict(dictionary)
        results_df.append(temp_df)
    '''


    #Substitutive code for single list taxid
    #--------------
    #Taxid with list
    blast_results = blastn_with_list(fasta_record.seq, query_size=query_size, list_taxid=taxid_list)
    #Parse handler to xml_string_to_handler to be able to feed it in SearchIO.read
    handler = xml_string_to_handler(blast_results)
    
    dictionary = blast_to_dictionary(handler)
    results_df = pd.DataFrame.from_dict(dictionary)


    print('Prepared DF')
    print(f'In {round(time.perf_counter()-start,2)}')
    print(len(results_df))
    #--------------




    #Filter DF by E-value / Bitscore / Identity
    filtered_blast_df = filter_df_blast(results_df, fasta_record, evalue=evalue_threshold, difference_from_query=len_threshold, identity_threshold=identity_threshold)
    print('Done filtering')
    print(f'In {round(time.perf_counter()-start,2)}')


    #Filter any futher? 


    #Retrieve all results from genbank
    #Retrieve accession list
    accession_list = filtered_blast_df['Accession'].tolist()
    #Retrieve genbank results, the handlers are managed in the function and return a JSON nested-type object
    genbank_results = retrieve_all_genbank(accession_list, 'A.N.Other@example.com')

    print('Done Genbank')
    print(f'In {round(time.perf_counter()-start,2)}')


    #Useless ATM
    '''    
    with open('genbank_results.txt', 'w') as savefile:
        savefile.write(str(genbank_results))
    '''

    #Create a dictionary with the genbank results
    genbank_dictionary = genbank_to_dictionary(genbank_results)

    #Create DF from genbank dictionary
    genbank_df = pd.DataFrame.from_dict(genbank_dictionary)

    print('Parsed Genbank')
    print(f'In {round(time.perf_counter()-start,2)}')


    #Merge DataFrames
    left = filtered_blast_df.loc[:,['Accession', 'ID', 'Description','Seq_length', 'Evalue', 'Bitscore', 'Tot_aln_span', 'Identity']]
    right = genbank_df.loc[:,['Accession', 'Accession_version','Organism', 'Taxonomy', 'Taxon', 'Nuc_sequence','Protein_ID', 'Prot_sequence']]
    combined_df = pd.merge(left, right, on='Accession')

    print('Merged dictionaries')
    print(f'In {round(time.perf_counter()-start,2)}')

    #------------
    #REMOVE
    #Save DF of sequences that are going to be aligned 
    with open(f'{output_df}_unfiltered.csv', 'w') as handle:
        handle.write(combined_df.to_csv())

    #-----------






    #Filter by taxon to prepare for alignment
    filtered_df = filter_df_taxon(combined_df, n_of_sequences=sequences_per_taxon)

    #Save DF of sequences that are going to be aligned 
    with open(output_df, 'w') as handle:
        handle.write(filtered_df.to_csv())


    #Declare protein
    protein_seq = None
    #Check if user specified ambiguities for AA translation
    if import_aa_sequence:
        protein_seq = open_fasta(query_protein)
    else:
        #Get protein from query if not specified otherwise
        protein_seq = fasta_record.translate(to_stop=True)

    


    #Prepare for multiple alingment
    fasta_string = fasta_for_alignment(protein_seq, filtered_df)


    #Run mafft alingment from string
    #mafft_alignment = run_mafft_from_string(fasta_string, mafft_directory)

    #Run mafft alignment saving file
    mafft_alignment = run_mafft_saving_file(fasta_string, mafft_directory, 'multiple_seq_fasta')

    print('Alignment done')
    print(f'In {round(time.perf_counter()-start,2)}')

    #Save mafft alignment
    with open(output_alignment, 'w') as savefile:
        savefile.write(mafft_alignment)
    
    
    
    
    #Open AlignIO from fasta file 

    alignment = AlignIO.read(output_alignment, 'fasta')

    #Construct tree from alignment #To change 
    tree = tree_from_alignment(alignment)

    #Write tree in xml
    Phylo.write(tree, output_xml_tree, 'phyloxml')

    #Draw tree
    with open('tree.ascii', 'w') as savefile:
        savefile.write(str(Phylo.draw_ascii(tree)))

    print('Done!')

    end = time.perf_counter()
    print(f'Finished in {round(end-start,2)}')



    '''

    Make fasta for alignment 

    '''












    '''
    #Write proteins for alignment in a file to feed to mofft
    with open('sequences_for_alignment_filtered_v2.fasta', 'w') as savefile:
        savefile.write(fasta_for_alignment)
    '''