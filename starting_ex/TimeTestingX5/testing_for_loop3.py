import schedule
import time
import subprocess
import os

print('imported packages')


def run_for_loop():
    scripts = ['multiprocessing_blast_v2_testing.py']

    print(f'Started {scripts}')
    script_list = ['/Users/Gioele/miniconda3/bin/python3', str(scripts)]      
    with open('testing_log.log', 'a') as handle:
        subprocess.run(script_list, stdout = handle)
    print(f'Finished running {scripts}')

            

print('scheduling')

schedule.every(1).hour.do(run_for_loop)

print('scheduled')

run_for_loop()

while not os.path.exists("/home/22/job_stop.txt"):
    schedule.run_pending()
    time.sleep(1)


print('terminated')