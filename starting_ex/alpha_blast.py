import subprocess

def run_blasts(input_file, taxdics, hls, thrds, b_res_filename):
    for t in list(taxdics):
        cmd = "/software/blast_v2.8.1/bin/blastp -db /home2/db/blast_v5/refseq_protein_v5 -query %s -taxids %s -max_target_seqs %s -outfmt 5 -num_threads %i -out %s" % (input_file, taxdics[t], hls, thrds, b_res_filename) 
        subprocess.Popen([cmd], shell = True,close_fds=True).wait()
