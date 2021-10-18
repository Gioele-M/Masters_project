!/Users/Gioele/miniconda3/bin/python3

import schedule
import time
import subprocess
import os


list_of_scripts = ['entrez_blast_testing.py']


def run_for_loop(list):
    for _ in range(5):    
        for script in list:        
            with open('testing_log.log', 'a') as handle:
                subprocess.run(script, stdout = handle)


schedule.every(1).minute.do(run_for_loop(list_of_scripts))

while not os.path.exists("/home/22/job_stop.txt"):
    schedule.run_pending()
    time.sleep(1)