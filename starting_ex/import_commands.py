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
command_parser.add_argument('-l', '--local_query', action='store_true', default=False) #action='store_true' -> returns true if flag is called 
command_parser.add_argument('-t', '--threading', action='store_true', default=False)

command_parser.add_argument('-e', '--evalue_threshold', type=float, default=10**-10)
command_parser.add_argument('-n', '--length_threshold', type=float, default=50.0)
command_parser.add_argument('-i', '--identity_threshold', type=float, default=50.0)

command_parser.add_argument('-s', '--sequences_per_taxon', type=int, default=1)
command_parser.add_argument('-s', '--sequences_per_taxon', type=int, default=1)


args = command_parser.parse_args()

print('{}'.format(args))
print('input: {}'.format(args.input_file)) 