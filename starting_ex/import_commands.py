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