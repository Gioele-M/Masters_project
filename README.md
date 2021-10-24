#BlastPTree
BlastPTree is a pipeline coded in Python 3.5.6 which provides a fast tool for preliminary homology analysis and gene tree inference. The pipeline uses external tools and packages including Biopython, Pandas, ETE3 and MAFFT and requires being operated in a conda environment. 


##Prerequisites:
The pipeline requires to have both *[MAFFT] (https://mafft.cbrc.jp/alignment/software/)* and *[conda] (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)* installed. Both softwares can be downloaded from the linked sources. 

##Installation:
The pipeline requires a conda environment to function properly. The setup YML file for the compatible environment is stored in the Git repository. The first step is cloning locally the repository.
‘git clone https://github.com/Username/final_branch.git'

Once the repository has been cloned, the YML file can be used to install the required conda environment.
‘’’sh
conda env create -f BlastPTree_environment.yml
‘’’
After installing the conda environment, the user should manually check if the installation succeeded typing
‘’’sh
conda BlastPTree_env list
‘’’

##Usage
The pipeline requires activating the conda environment
‘’’sh
‘conda activate BlastPTree_env’
‘’’

The pipeline has four compulsory positional arguments and six optional arguments
Compulsory arguments
i) the input string (input string can either be a nucleotide/protein fasta file or an accession number)
ii) the MAFFT directory 
iii) the user’s email address (required by the biopython Bio.Entrez.efetch() module to connect to the NCBI servers) 
iv) the output name

Optional arguments
‘-e/--evalue_threshold’ sets the threshold for the E-value for Blast results (default 10^-10); 
‘-n/--length_threshold’ sets the threshold for the length of the retrieved sequences. The tolerance is based on the query length at which it is added and subtracted X% (default 50). For example, if the query is 100 sequences long, only results long between 50 and 150 are considered.
‘-i/--identity_threshold’ sets the cutoff value for the minimum allowed sequence identity (default 50%). 
‘-x/--taxid_list’ inputs the requested list of taxid (default None); 
-b/--retrieve_blast sets the maximum number of results downloaded from Blast (default 100); 
'-s/--sequences_per_taxon’ sets the maximum number of sequences per taxon included in the final output (default 1). 


To use the pipeline, the Python interpreter should be called first, followed by the BlastPTree.py script and its arguments. A standard use of the pipeline is:

‘’’sh
python BlastPTree.py input_string /MAFFT/directory/mafft example@email.com output -e 10^-10 -n 50 -i 50 -x 9595 9564 -s 1 -b 100
‘’’



##Contact
Gioele Montis - gioelemontis97@gmail.com
