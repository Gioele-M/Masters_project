import argparse

#USER INPUTS
filename = 'sry_gene.fasta' #'human_mx1.fas'
mafft_directory = r'/Users/Gioele/miniconda3/bin/mafft'
email = 'A.N.Other@example.com'
output_name = 'run1'

taxid_list = [] #['9592']

#OPTIONAL USER INPUTS
local_query = False
threading = False
query_size = 100
server, database = '',''#FOR BLAST

evalue_threshold = 10**-10
len_threshold = 50
identity_threshold = 50

sequences_per_taxon = 3
import_aa_sequence = False
query_protein = 'sry_protein.fasta'


type_of_molecule = 'nucleotide'





command_parser = argparse.ArgumentParser(description='Run your BLAST and obtain MSA and phylogenetic tree')

command_parser.add_argument('input_file')
command_parser.add_argument('mafft_directory')
command_parser.add_argument('email')
command_parser.add_argument('output_name')

#Taxid list as optional list
command_parser.add_argument('-x', '--taxid_list', nargs='*', type=str)

#optional argument
command_parser.add_argument('-e', '--evalue_threshold', type=float, default=10**-10)
command_parser.add_argument('-n', '--length_threshold', type=float, default=50.0)
command_parser.add_argument('-i', '--identity_threshold', type=float, default=50.0)

command_parser.add_argument('-q', '--blast_query_size', type=int, default=100)
command_parser.add_argument('-s', '--sequences_per_taxon', type=int, default=1)


args = command_parser.parse_args()

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


print(input_string, mafft_directory, email, output_name, taxid_list, evalue_threshold, len_threshold, identity_threshold, query_size,sequences_per_taxon)
