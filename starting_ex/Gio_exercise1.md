# Starting with python & BLAST

Exercises to get you started with using python and biopython + become aware of sequence identity search algorithms 


## Learn the python basics (lists, dictionaries, etc)
- [basics python tutorial](https://www.learnpython.org/en/Welcome)
- [play around with pandas](https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html)
- [get started with biopython](https://www.tutorialspoint.com/biopython/biopython_introduction.htm)


## Write a BLAST search python script

- import coding sequence from fasta file with biopython's [SeqIO.parse](https://biopython.org/wiki/SeqIO)
- translate the nucleotide sequence to amino acid with [Bio.Seq](https://biopython.org/docs/1.75/api/Bio.Seq.html)
- perform both blastn and blastp with the respective sequences using the [NCBI biopython module](https://biopython.org/docs/1.75/api/Bio.Blast.NCBIWWW.html)
	- can you figure out how to perform the BLAST search for specific species by [taxid](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)? 
- save the blast output  results as an xml and parse the results back to python with [Bio.SearchIO](https://biopython.org/docs/1.75/api/Bio.SearchIO.html)
- once you've played around with what's in the blast output and how to parse it try to summarise it in a table using pandas:
	- [turn a dictionary into a dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.from_dict.html)
	- [sort by any column you want](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sort_values.html)
	- [save to a csv file (without the index!)](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html)
- can you use [input](https://www.w3schools.com/python/ref_func_input.asp) to run your script on any fasta file?



## Compile a list of sequence similarity methods (like BLAST)
- have a look around google for sequence similarity methods similar to blast that can be used in python
- are there any interesting alternatives?
- maybe some machine learning based approaches?
- is it easy to access available databases (e.g. Genbank or Ensembl) with these methods?
