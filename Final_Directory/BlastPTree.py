import tempfile
import subprocess
import logging
import time
import traceback
import argparse
import pandas as pd
from pandas.core.frame import DataFrame
from Bio import SeqIO, SeqRecord, SearchIO, AlignIO, Phylo
from Bio.Blast import NCBIWWW
import Bio.Entrez
from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from ete3 import PhyloTree


#Function to open fasta file, returns a SeqRecord object
def open_fasta(filename) -> SeqRecord:
    with open(filename) as handle:
        sequence_record = SeqIO.read(handle, 'fasta')
    logging.info('Opened sequence {}'.format(sequence_record.id))
    return sequence_record


#Function to run BLAST with taxid list
def blastp_with_list(sequence, list_taxid = [], query_size = 200):
    result_handler, result_storer = None, None
    #If list is empty run query without specific taxid
    if len(list_taxid) <1:
        result_handler = NCBIWWW.qblast('blastp', 'nr', sequence, hitlist_size=query_size) 
        result_storer = result_handler.read()
    #Otherwise prepare string of Entrez and parse it to qblast
    else:
        entrez_query = ''
        for taxid in list_taxid:
            entrez_query += 'txid{}[ORGN]'.format(taxid)
            if taxid != list_taxid[-1]:
                entrez_query += ' OR '
        result_handler = NCBIWWW.qblast('blastp', 'nr', sequence, entrez_query= entrez_query, hitlist_size=query_size)
        result_storer = result_handler.read()
    logging.info('BLASTp specifying {} taxid(s) completed'.format(len(list_taxid)))
    #Returns Blast XML string
    return result_storer


#Function to create a temporary file from string and read it with SearchIO.read
def xml_string_to_handler(string):
    #make temporary file
    tmp = tempfile.NamedTemporaryFile(mode='a+')
    #write string
    tmp.write(string)
    #Read blast results
    handler = SearchIO.read(tmp.name, 'blast-xml')
    #Close temporary file
    tmp.close()
    logging.info('Blast returned {} results'.format(len(handler)))
    return handler


#Creation of a dictionary from the Blast results
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
        #Collect data of all hsps for each hit
        for hsp in result.hsps:
            all_alnspan.append(int(hsp.aln_span))
            all_gapnum.append(int(hsp.gap_num))
        #Calculate total alignment span and gaps to calculate identity
        tot_alnspan, tot_gapnum = int(), int()
        seq_len = int(result.seq_len)
        for span in all_alnspan:
            tot_alnspan += span
        for gap in all_gapnum:
            tot_gapnum += gap
        identity = ((tot_alnspan - tot_gapnum)/seq_len)*100
        blast_dictionary['Tot_aln_span'].append(tot_alnspan)
        blast_dictionary['Identity'].append(round(identity, 3))

    logging.info("{} entries were recorded from the BLASTp results".format(len(blast_dictionary['ID'])))
    #Returns dictionary
    return blast_dictionary


#Filter DF based on Evalue, sequence length and identity
def filter_df_blast(df:DataFrame, query:SeqRecord, evalue = 10**-10, difference_from_query = 50, identity_threshold=50) -> DataFrame:
    #Record initial dataframe length
    initial_len = len(df)

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
    percent_query = query_length*(difference_from_query/100)
    lower, upper = query_length-percent_query, query_length+percent_query
    #Filter DF
    df_to_return = df_to_return[df_to_return['Seq_length'] > lower]
    df_to_return = df_to_return.reset_index(drop=True)
    df_to_return = df_to_return[df_to_return['Seq_length'] < upper]
    df_to_return = df_to_return.reset_index(drop=True)

    #Filter by identity
    df_to_return = df_to_return[df_to_return['Identity'] > identity_threshold]
    df_to_return = df_to_return.reset_index(drop=True)
    
    #Calculate how many sequences were removed
    final_len = len(df_to_return)

    logging.info("{} sequences were removed filtering the BLAST results DF, returning a DF with {} entries".format((initial_len-final_len),final_len))
    return df_to_return


#Retrieve all results' entries from efetch
def retrieve_all_efetch(list_of_entries, email):
    Bio.Entrez.email = email
    list_of_sequences = []
    for entry in list_of_entries:
        handler = Bio.Entrez.efetch(db='protein', id=entry, rettype = 'fasta',retmode = 'xml', retmax=1) #Returns JSON regardless
        gb_info = Bio.Entrez.read(handler, 'text')#Returns nested lists and dictionaries 
        list_of_sequences.append(gb_info)
    logging.info('{} protein entries were retrieved from NCBI protein database'.format(len(list_of_sequences)))
    #Returns list of JSON strings
    return list_of_sequences


#Efetch to dictionary parser
def efetch_protein_to_dictionary(list_of_efetch):
    #Declare new dictionary
    dictionary = {'Accession':[],'Protein_ID':[], 'Taxid':[], 'Organism_name':[], 'Description':[], 'Seq_length':[], 'Prot_sequence':[]}
    #Determine if all parsing was successful
    error_message = False
    not_parsed = 0
    for wrapper in list_of_efetch:
        try:
            #Cast into dictionary to avoid random exceptions
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
        except KeyError:
            #Flag if an exception is raised so to be reported in the log file
            error_message = True
            not_parsed += 1        
    if error_message:
        logging.error('{} sequence(s) failed parsing from EFetch'.format(not_parsed))
    logging.info("{} sequences were parsed correctly from EFetch".format(len(dictionary['Accession'])))
    return dictionary


#Function to filter the final DF by taxon
def filter_df_taxon(df:DataFrame, n_of_sequences = 1) -> DataFrame:
    #Record initial length of DF
    initial_len = len(df)
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
        #Sort it by value and only get the X first sequences
        temp_df = temp_df.sort_values(['Evalue', 'Identity', 'Bitscore'], ascending=[True, False, False]) #('Identity', ascending=False)
        if len(temp_df) > 0:
            if len(temp_df) > n_of_sequences:
                df_toreturn = df_toreturn.append(temp_df[:n_of_sequences])
            else:
                df_toreturn = df_toreturn.append(temp_df)

    #Refactor indexes
    df_toreturn = df_toreturn.reset_index(drop=True)
    #Calculate final length and append info to logging 
    final_len = len(df_toreturn)
    logging.info("Filtering by taxon: {} unique taxid(s) collected, {} entries discarded".format(len(retrieved_taxids), (initial_len - final_len)))
    return df_toreturn


#Prepare string for alingment file
def fasta_for_alignment(query:SeqRecord, df:DataFrame) -> str:
    #Initialise string and add query
    string = ''
    string += ">Query:{}\n".format(query.id)
    string += "{}\n".format(query.seq)
    #Add all sequence IDs and amino-acid sequences
    for i in range(len(df['Accession'])):
        #Check if protein is missing
        if len(df['Prot_sequence'][i]) <1:
            logging.error('Sequence {} could not be added to the MSA because missing the sequence'.format(df['Accession'][i]))
            continue
        string += ">{}\n".format(df['Protein_ID'][i])
        string += "{}".format(df['Prot_sequence'][i])
        #Check if is not the last element in the list, if not add a return
        if i != len(df['Accession']):
            string += '\n'
    logging.info("{} sequences were prepared for alignment with the query".format(len(df)))
    return string


#Run mafft by saving the fasta sequences for alignment in a file and passing it to mafft 
def run_mafft_saving_file(fasta:str, mafft_directory:str, filename:str) -> str:
    #Write fasta file for alignment
    file = '{}.fasta'.format(filename)
    with open(file, 'w') as handle:
        handle.write(fasta)
    #Parse and run with mafft 
    command_list = [mafft_directory, '{}'.format(file)] 
    process = subprocess.Popen(command_list, universal_newlines= True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()
    #Save standard error of mafft
    with open('{}_mafft.stderr'.format(filename), 'w') as handle:
        handle.write(stderr)
    #Return alignment 
    return stdout


#Function to get protein sequence and save fasta file given an accession number
def get_fasta_from_accession(accession, email):
    #Get the efetch handler for request accession
    Bio.Entrez.email = email
    handler = Bio.Entrez.efetch(db='protein', id=accession, rettype = 'fasta',retmode = 'xml', retmax=1) 
    query_protein_efetch = Bio.Entrez.read(handler, 'text')
    
    #Make a dictionary -- passed as list to recycle the efetch_protein_to_dictionary function
    dictionary_query = efetch_protein_to_dictionary([query_protein_efetch])

    #Make and save fasta file of the query
    fasta_string_query = ">{} \n{}".format(dictionary_query['Accession'][0], dictionary_query['Prot_sequence'][0])
    fasta_file_name = "{}_sequence.fasta".format(dictionary_query['Accession'][0])
    with open(fasta_file_name, 'w') as handle:
        handle.write(fasta_string_query)

    #Open sequence with open_fasta and return it 
    fasta_record_q = open_fasta(fasta_file_name)
    return fasta_record_q


#Function to open input string
def open_input(input, email):
    #Declare empty protein sequence
    protein_sequence = None
    #Check that the parsed string has a fasta extension, if it does, pass it to open fasta function
    if input.endswith('.fas') or input.endswith('.fasta'):
        sequence = open_fasta(input)
        #Try translating the sequence, if an error is raised the sequence is already a peptide and can be returned
        try: 
            protein_sequence = sequence.translate(to_stop = True)
            logging.info('Nucleotide sequence was translated to protein')
        except Exception:
            protein_sequence = sequence
            logging.info('Protein sequence was opened')
    #If the string doesn't have an extension it will be passed to the get_fasta_from_accession function which produces a fasta file
    else:
        try:
            protein_sequence = get_fasta_from_accession(input, email)
            logging.info('Protein sequence was retrieved from NCBI protein database')
        except Exception as e:
            logging.error(traceback.format_exc())
            logging.error('The file must have a .fas or .fasta extension')

    return protein_sequence


#Function to create a gene tree from the MSA
def tree_from_alignment(alignment):
    #Create a calculator for identity and distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    #Build the tree with the neighbour joining method
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    logging.info('Tree was produced using the neighbour joining method')
    return tree


if __name__ == '__main__':
    #Import user inputs
    command_parser = argparse.ArgumentParser(description='Run your BLAST and obtain MSA and phylogenetic tree')
    #Compulsory arguments
    command_parser.add_argument('input_file')
    command_parser.add_argument('mafft_directory')
    command_parser.add_argument('email')
    command_parser.add_argument('output_name')
    #Taxid list as optional list
    command_parser.add_argument('-x', '--taxid_list', nargs='*', type=str)
    #Optional arguments
    command_parser.add_argument('-e', '--evalue_threshold', type=float, default=10**-10)
    command_parser.add_argument('-n', '--length_threshold', type=float, default=50.0)
    command_parser.add_argument('-i', '--identity_threshold', type=float, default=50.0)
    command_parser.add_argument('-q', '--blast_query_size', type=int, default=100)
    command_parser.add_argument('-s', '--sequences_per_taxon', type=int, default=1)
    #Get the arguments
    args = command_parser.parse_args()
    #Assign them to the respective variables
    input_string = args.input_file
    mafft_directory = args.mafft_directory #r'/Users/Gioele/miniconda3/bin/mafft'
    email = args.email
    output_name = args.output_name
    taxid_list = [] #['9592']
    if type(args.taxid_list) == list:
        taxid_list = args.taxid_list
    evalue_threshold = args.evalue_threshold
    len_threshold = args.length_threshold
    identity_threshold = args.identity_threshold
    query_size = args.blast_query_size
    sequences_per_taxon = args.sequences_per_taxon


    #Declare logging file
    logging.basicConfig(filename='{}.log'.format(output_name), filemode='w', format='%(levelname)s:%(message)s', level=logging.DEBUG) #logging refreshes every run and only displays type of message: message

    #Filenames for outputs
    output_df = '{}_df.csv'.format(output_name)
    output_alignment = '{}_alignment.fasta'.format(output_name)
    output_xml_tree = '{}_xml_tree.xml'.format(output_name)
    output_tree_newick = '{}_newick_tree.nwk'.format(output_name) #!!!!
    output_tree_png = '{}_tree_image.png'.format(output_name)


    #Start of pipeline
    start = time.perf_counter()

    #Open fasta record
    fasta_record = open_input(input_string, email)    
    print('Opened input record in {} seconds'.format(round(time.perf_counter()-start, 2)))

    #Run BLAST
    blast_results = blastp_with_list(fasta_record.seq, query_size=query_size, list_taxid=taxid_list)
    #Parse BLAST results to temporary file so SearchIO.read can be used
    handler_blast = xml_string_to_handler(blast_results)
    #Make BLAST dictionary from handler
    dictionary_blast = blast_to_dictionary(handler_blast)
    #Cast dictionary into DF
    results_df_blast = pd.DataFrame.from_dict(dictionary_blast)
    print('Prepared BLAST DF with {} entries in {} seconds'.format(len(results_df_blast), round(time.perf_counter()-start, 2)))
    #Filter DF by E-value / Bitscore / Identity
    filtered_blast_df = filter_df_blast(results_df_blast, fasta_record, evalue=evalue_threshold, difference_from_query=len_threshold, identity_threshold=identity_threshold)
    print("BLAST DF after filtering holds {} entries".format(len(filtered_blast_df)))


    #Retrieve all BLAST results' sequences and information with efetch()
    #Retrieve protein accession list
    accession_list = filtered_blast_df['Accession'].tolist()
    #Get information on sequences from the NCBI servers using the efetch() function
    efetch_result = retrieve_all_efetch(accession_list, email)
    # Make dictionary from efetch() results
    dictionary_efetch = efetch_protein_to_dictionary(efetch_result)
    #Cast dictionary into DF
    efetch_df = pd.DataFrame.from_dict(dictionary_efetch)
    print('Efetch retrieved and parsed {} entries in {} seconds'.format(len(efetch_df), round(time.perf_counter()-start,2)))


    #Merge BLAST and efetch DataFrames
    left = filtered_blast_df.loc[:,['Accession', 'ID','Seq_length', 'Evalue', 'Bitscore', 'Tot_aln_span', 'Identity']]
    right = efetch_df.loc[:,['Accession', 'Protein_ID','Taxid', 'Organism_name', 'Description', 'Prot_sequence']]
    combined_df = pd.merge(left, right, on='Accession')
    #Filter DF based on taxons
    filtered_df = filter_df_taxon(combined_df, n_of_sequences=sequences_per_taxon)
    print("The combined DF having {} entries was reduced to {} entries after keeping {} sequence per taxid".format(len(combined_df), len(filtered_df), sequences_per_taxon))
    print('Filtered DF presents {} entries'.format(len(filtered_df)))


    #Output DF of sequences that are going to be aligned in CSV format
    filtered_df.to_csv(output_df, index = False)
    print('DataFrame was output in CSV format in file {}'.format(output_df))
    #Output TSV for FigTree
    tsv_df = filtered_df[['Protein_ID'] + [col for col in filtered_df.columns if col!= 'Protein_ID']]
    tsv_df = tsv_df.drop('Prot_sequence', axis = 1)
    tsv_df.to_csv('{}.tsv'.format(output_df), sep = '\t', index = False)
    print('DataFrame compatible with FigTree was output in file {}.tsv'.format(output_df))


    #Prepare fasta file for multiple sequence alignmetn
    fasta_string_for_alignment = fasta_for_alignment(fasta_record, filtered_df)
    #Run mafft alignment
    mafft_alignment = run_mafft_saving_file(fasta_string_for_alignment, mafft_directory, output_name)
    print('Alignment with mafft was produced in {} seconds'.format(round(time.perf_counter()-start,2)))
    #Save mafft alignment
    with open(output_alignment, 'w') as savefile:
        savefile.write(mafft_alignment)
    print('Mafft alignment was saved in file {}'.format(output_alignment))


    #Open alignment from fasta file using AlignIO
    alignment = AlignIO.read(output_alignment, 'fasta')
    #Construct tree from alignment 
    tree = tree_from_alignment(alignment)  

    #Output tree in PhyloXML format
    Phylo.write(tree, output_xml_tree, 'phyloxml')
    #Output tree in Newick format
    Phylo.write(tree, output_tree_newick, 'newick')
    print('The inferred tree was output in PhyloXML and Newick format in files {} and {}'.format(output_xml_tree, output_tree_newick))


    #Re-import newick string to feed to ETE3
    with open(output_tree_newick, 'r') as handle:
        newick_string = handle.read()

    #Create a DF with only protein ids and organisms names and drop the index
    #This is required to change the leaf nodes in the ETE3 final output image
    dict_of_df = filtered_df[['Protein_ID', 'Organism_name']]
    dict_of_df.reset_index(drop=True, inplace=True)
    #Create a dictionary from DF 
    dict_of_df = dict_of_df.to_dict('records')
    
    #Replace all organisms names' spaces with underscores to avoid confusion when using the names
    for entry in dict_of_df:
        entry['Organism_name'] = str(entry['Organism_name']).replace(' ', '_')

    #Change all instances of protein ID to Organism name | protein ID in both newick and MSA strings
    #Makes the ETE3 output more human readable
    for entry in dict_of_df:
        #Get new strings to be applied
        pro_id = str(entry['Protein_ID'])
        org_name = str(entry['Organism_name'])
        new_entry = '{}|{}'.format(org_name, pro_id)
        #Change strings in both Newick tree and MSA results
        newick_string= newick_string.replace(pro_id, new_entry)
        mafft_alignment= mafft_alignment.replace(pro_id, new_entry)


    #Create PhyloTree object from newick string
    msa_tree = PhyloTree(newick_string, quoted_node_names=True, format=1) 

    #Rooting of the tree
    #Returns the node that divides the current tree into two distance-balanced partitions.
    R = msa_tree.get_midpoint_outgroup()
    #Sets a descendant node as the outgroup of a tree
    msa_tree.set_outgroup(R)

    #Link tree to the MAFFT Alignment 
    msa_tree.link_to_alignment(alignment=mafft_alignment, alg_format='fasta')
    #Render final tree
    msa_tree.render(output_tree_png)
    print('Inferred tree image was saved in file {}'.format(output_tree_png))

    end = time.perf_counter()
    print('Run has completed in {} seconds'.format(round(end-start,2)))
    print('Inferred tree was produced comprehending {} sequences. \nCheck the log file {}.log for detailed results'.format(len(filtered_df), output_name)) 