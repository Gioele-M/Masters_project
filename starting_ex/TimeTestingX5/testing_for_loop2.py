import schedule
import time
import subprocess
import os

print('imported packages')


def run_for_loop():
    list_of_scripts = ['multithreading_blast_v1.1_testing.py']

    for x in range(5):    
        for script in list_of_scripts:  
            print(f'Started {script}')
            script_list = ['/Users/Gioele/miniconda3/bin/python3', str(script)]      
            with open('testing_log.log', 'a') as handle:
                subprocess.run(script_list, stdout = handle)
            print(f'Finished running {x+1} {script}')
            time.sleep(5)

            

print('scheduling')

schedule.every(1).hour.do(run_for_loop)

print('scheduled')

run_for_loop()

while not os.path.exists("/home/22/job_stop.txt"):
    schedule.run_pending()
    time.sleep(1)


print('terminated')