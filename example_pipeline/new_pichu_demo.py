#### You can run this yourself on alpha! (I think...) ####
#### just make sure to create an empty folder named:  ####
#### 'run' + the run_id in the ctl file, e.g. run001  ####
#### Once fasta, taxa txt, ctl + run_id subfolder are ####
#### in the same directory, run as:                   ####
        #### nohup python3 new_pichu_demo.py &        ####



#import dependencies
from Bio.Blast import NCBIWWW
from Bio import SearchIO
from Bio import SeqIO
import pandas as pd
from Bio import Entrez
import subprocess
import xml.etree.ElementTree as ET
import os



def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]


#### This is made to run on alpha with a local blast db ####

def run_blasts(taxdics, theref):
    for t in list(taxdics):
        cmd = "/software/blast_v2.8.1/bin/blastp -db /home2/db/blast_v5/refseq_protein -query %s -taxids %s -max_target_seqs 1000 -outfmt 5 -out %s.xml" % (taxdics[t], theref, t) 
        subprocess.Popen([cmd], shell = True,close_fds=True).wait()



#for parsing the ctl xml file info#
def xmltext(root, var):
    for child in root: 
        if child.tag == var:
            return child.text

def xmllist(root, var):
    for child in root: 
        if child.tag == var:
            thelist = []
            for schild in child:
                if schild.text == 'True':
                    thelist.append(schild.tag)
            return thelist
####################################






tree = ET.parse('ctl.xml')

root = tree.getroot()


dir = xmltext(root, 'dir')

organisms = xmltext(root, 'organisms')




####    CONTROL    ####

#NCBI requires email
Entrez.email = xmltext(root, 'email')

# print(Entrez.email)

#to keep track of different runs
runnum = xmltext(root, 'run_id')

# print(runnum)

#this should be defined in the control
#use a specific reference for everything - path relative to py script
#will use mkdir from here
input_file = xmltext(root, 'input')

# print(input_file)

#define the number of hits outputed from the BLAST
hls = int(xmltext(root, 'b_hits'))

# print(hls)


#define number of threads for blast
thrds = int(xmltext(root, 'threads'))


#define the e value filter cutoff
e_val_filter = float(xmltext(root, 'e_val_filter'))

# print(e_val_filter)

#TD# 
#have an option to select which outputs to get (default all)
theoutputs = xmllist(root, 'outputs')

# print(theoutputs)


########################


#I will append all the final species dfs here
alldflist = []


#create a dictionary with species names and corresponding taxids based on the ctl file#
taxdics = {}

with open(organisms, 'r+') as thetaxa:
    for line in thetaxa:
        line = line.replace('\n', '')
        orgname = line.replace(' ', '-')
        thetaxid = get_tax_id(line)
        taxdics.update({orgname:thetaxid})


#export taxdics to use for future runs
with open('run%s_taxids.txt'%runnum, 'w+') as out:
    out.write(str(taxdics))


# #cut down on time by having this premade instead of running the above block
    #the Haplorrhini version is also available in a text file in the example_pipeline folder
# taxdics = {'Homo-sapiens': '9606', 'Gorilla-gorilla-gorilla': '9595', 'Macaca-mulatta': '9544', 'Pan-troglodytes': '9598', 'Cercopithecus-mona': '36226', 'Alouatta-palliata': '30589', 'Cercopithecus-neglectus': '36227', 'Pithecia-pithecia': '43777', 'Saguinus-imperator': '9491', 'Ateles-geoffroyi': '9509', 'Piliocolobus-tephrosceles': '591936', 'Carlito-syrichta': '1868482', 'Mandrillus-leucophaeus': '9568', 'Trachypithecus-francoisi': '54180', 'Erythrocebus-patas': '9538', 'Mandrillus-sphinx': '9561', 'Pygathrix-nemaeus': '54133', 'Theropithecus-gelada': '9565', 'Plecturocebus-donacophilus': '230833', 'Aotus-nancymaae': '37293', 'Cercocebus-atys': '9531', 'Macaca-nemestrina': '9545', 'Chlorocebus-sabaeus': '60711', 'Hylobates-moloch': '81572', 'Sapajus-apella': '9515', 'Pan-paniscus': '9597', 'Pongo-pygmaeus': '9600', 'Rhinopithecus-bieti': '61621', 'Semnopithecus-entellus': '88029', 'Rhinopithecus-roxellana': '61622', 'Nasalis-larvatus': '43780', 'Cebus-albifrons': '9514', 'Saimiri-boliviensis-boliviensis': '39432', 'Macaca-fascicularis': '9541', 'Nomascus-leucogenys': '61853', 'Callithrix-jacchus': '9483', 'Papio-anubis': '9555', 'Pongo-abelii': '9601'}

########################



#error file text
errorstr = 'Cannot find in BLASTdb:'


#add the loop
for t in list(taxdics):

    animal = t


    ####    filenames    ####

    #raw BLAST results xml
    # b_res_filename = str(animal + '/my_blast' + runnum + '.xml')
    b_res_filename = 'run%s/%s.xml'%(runnum, animal)
    
    #unfiltered BLAST results table
    # bdata_all_out = str(animal + '/' + animal + '_blast' + runnum + '_meta_all.tsv')
    bdata_all_out = 'run%s/%s_blast%s_meta_all.tsv'%(runnum, animal, runnum)

    #filtered BLAST results table
    # bdata_out = str(animal + '/' + animal + '_blast' + runnum + '_meta.tsv')
    bdata_out = 'run%s/%s_blast%s_meta.tsv'%(runnum, animal, runnum)

    #hits protein fasta file
    # outfnprot = str(animal + '/' + animal + '_BTNs_' + runnum + '_prot.fas')
    outfnprot = 'run%s/%s_hits%s_prot.fas'%(runnum, animal, runnum)

    #hits protein genprot file
    # outfnprot = str(animal + '/' + animal + '_BTNs_' + runnum + '_prot.gp')
    outfnprotgp = 'run%s/%s_hits%s_prot.gp'%(runnum, animal, runnum)

    #hits coding sequence fasta file
    # outfnnt = str(animal + '/' + animal + '_BTNs_' + runnum + '_cds.fas')
    outfnnt = 'run%s/%s_hits%s_cds.fas'%(runnum, animal, runnum)

    ########################

    
    
    #read reference
    read_seq = SeqIO.read(input_file, 'fasta')

    #query length, the length of the reference
    qlen = len(read_seq)

    #store as fasta string for qblast
    seq_dat = str('>' + read_seq.description + '\n' + str(read_seq.seq))    
    
    
    ####    alpha BLAST    ####
    
    #run the blast as an external bash command
    cmd = "/software/blast_v2.8.1/bin/blastp -db /home2/db/blast_v5/refseq_protein_v5 -query %s -taxids %s -max_target_seqs %s -outfmt 5 -num_threads %i -out %s" % (input_file, taxdics[t], hls, thrds, b_res_filename) 
    # print('%s\n\n'%cmd)
    theblast = subprocess.Popen([cmd],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True,close_fds=True)
    theblast.wait()
    if theblast.returncode != 0:
        # print("\n\n!! Can't find %s - %s in blastdb\n\n"%(t, taxdics[t]))
        errorstr = errorstr + '\n%s\t%s'%(t, taxdics[t])
    
    
    ###########################
    
    
    # #from here with ready xml
    
    
    #first check if the xml gave any results
    if os.stat(b_res_filename).st_size != 0:


        # #redefine b_res_filename:
        # b_res_filename = 'alltaxa.xml'
        #parse the hits from the BLAST xml
        blast_qresult = SearchIO.read(b_res_filename, "blast-xml")



        #create the dictionary with BLAST results
        blast_resdi = {'name':[], 'id':[], 'q_cover':[], 'e_val':[], 'hit_length':[]}

        for hit in blast_qresult:
                
            #the hit information
            hitdit = hit[0]
            
            #the hit id cleaned up
                #this might need further cleaning depending on the case
            hitid = str(hitdit.hit_id).replace('ref|', '')
            hitid = hitid.replace('|', '')
            
            #this q_cover metric is: 
                #the length of matching sequence between reference and hit over the length of the reference
            query_cov = (hitdit.query_span/qlen)

            #the length of the matching sequence
            hitlen = len(hitdit.hit) 

            
            #the cutoff here is unnecessary since it can be applied in a later stage
                #if query_cov > 0.5:
              
            
            #add everything to the dict
            
            blast_resdi['name'].append(str(hitdit.hit_description))

            blast_resdi['id'].append(hitid)

            blast_resdi['q_cover'].append(query_cov)
            
            blast_resdi['e_val'].append(hitdit.evalue)

            blast_resdi['hit_length'].append(hitlen)
            


        #create a pandas df
        df = pd.DataFrame.from_dict(blast_resdi)


        #save df to csv
        df.to_csv(bdata_all_out, sep = '\t', index = False)
        
        # print(df)




#### You don't need to read the below in detail, ####
#### but it might give you some ideas about info ####
####           to parse from the xml             ####

        ####    Remove Duplicates    ####

        df = pd.read_csv(bdata_all_out, sep = '\t')
        
        
        #only do for df that are not empty
        if len(df) != 0:

            #remove 'low quality protein' and 'predicted' tags
            df['name'] = df.name.str.replace('LOW QUALITY PROTEIN: ', '')
            df['name'] = df.name.str.replace('PREDICTED: ', '')
            
            
            #split the isoform bit (only works if the word isoform is in at least one of the rows)
            if len([x for x in list(df.name) if ' isoform' in x]) > 0:
                df[['name','isoform']] = df.name.str.split(" isoform",n=1, expand=True)
            else:
                df['isoform'] = ''


            #split the precursor bit (only works if the word isoform is in at least one of the rows)
            if len([x for x in list(df.name) if ' precursor' in x]) > 0:
                df[['name','precursor']] = df.name.str.split(" precursor",n=1, expand=True)
            else:
                df['precursor'] = ''


            #split the extra name bit (only works if the word isoform is in at least one of the rows)
            if len([x for x in list(df.name) if ' [' in x]) > 0:
                print([x for x in list(df.name) if ' [' in x])
                df[['name','extra']] = df.name.str.split(" \[", n=1, expand=True)
            else:
                df['extra'] = ''

            df.sort_values(by=['hit_length', 'isoform'], ascending = [False, True], inplace = True)


            df.drop_duplicates(subset='name', keep="first", inplace = True)


            #filter by e value
            df = df[df['e_val']<= e_val_filter]

            df.sort_values(by=['q_cover'], ascending = False, inplace = True)



            #export the df filtered by e value and for duplicates, sorted by q value
            df.to_csv(bdata_out, sep = '\t', index = False)

            #store all the protein ids
            idez = list(df.id)
            
            
            print('\n')
            print(df)
            print('\n\n')
            
            
            
            #write df with everything
            df['species'] = animal
            df['taxid'] = taxdics[t]
            df.drop(columns=['precursor' ,'extra'])
            alldflist.append(df)
            
            

#### There's 4 output files (depending on the options)      ####
#### 1. protein fasta for each hit                          ####
#### 2. coding sequence fasta for each hit                  ####
#### 3. protein genpept file for each hit - tad unnecessary ####
#### 4. final results csv file                              ####

            ####    PROTEIN FASTA    ####
            
            if 'protfas' in theoutputs:

                #efetch fetches the genome information in fasta format
                net_handle = Entrez.efetch(db="protein", id=idez, rettype="fasta", retmode="text")

                #opens a writable file with the output file name
                with open(outfnprot, "w") as out_handle:
                    out_handle.write(net_handle.read())
                    out_handle.close()
                    net_handle.close()



            ####    PROTEIN GENPEPT    ####
            
            if 'protgp' in theoutputs:

                #efetch fetches the genome information in fasta format
                net_handle = Entrez.efetch(db="protein", id=idez, rettype="gp", retmode="text")

                #opens a writable file with the output file name
                with open(outfnprotgp, "w") as out_handle:
                    out_handle.write(net_handle.read())
                    out_handle.close()
                    net_handle.close()


            ####    CDS FASTA    ####
            
            if 'cdsfas' in theoutputs:

                #efetch fetches the genome information in fasta format
                net_handle = Entrez.efetch(db="protein", id=idez, rettype="fasta_cds_na", retmode="text")

                #opens a writable file with the output file name
                with open(outfnnt, "w") as out_handle:
                    out_handle.write(net_handle.read())
                    out_handle.close()
                    net_handle.close()



            


#write error output
with open('errors.out', 'w+') as out:
    out.write(errorstr)


#write final all species df
finaldf = pd.concat(alldflist)
finaldf.to_csv('blast%s_fn.csv'%runnum, index = False)
print(finaldf)

