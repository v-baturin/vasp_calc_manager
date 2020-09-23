#!python3
import sys
import argparse
import socket
import os
from datetime import datetime
import pickle
import time
from src.common_tools import cjson_load
from src.tasks_database import GauCalcDB

# sys.stdout = open('log', 'w')

# if len(sys.argv) == 1:
#     sys.argv.extend(['-p', '1-10_even_POSCARS',
#                      '-f', 'RECALC',
#                      '-d', 'database.pkl',
#                      '-c', 'local', '--maxcalc', '3',
#                      '--min_mult'])
print(sys.argv)  # debug

ap = argparse.ArgumentParser()
ap.add_argument("-c", "--machine", required=False, help='Which cluster')
ap.add_argument("-f", "--recalc_folder", required=False, help='Calc folder', default='RECALC')
ap.add_argument("-p", "--poscars_file", required=False, help="POSCARS-file", default='POSCARS')
ap.add_argument("-d", "--database_file", required=False, help='DB file', default='database.pkl')
ap.add_argument("-m", "--maxcalcs", required=False, help='Maximum parallel calcs', default='1')
ap.add_argument("-o", "--outfile_pattern", required=False, help='Maximum parallel calcs', default='log')
ap.add_argument("-s", "--sleep", required=False, help='Sleep time, in minutes', default='5')
ap.add_argument("--machines", required=False, help='Machines cjson file', default='machines.cjson')
ap.add_argument("--scenarios", required=False, help='Scenarios cjson file', default='scenarios.cjson')
ap.add_argument("--min_mult", required=False, help='Use minimal multiplicity', action='store_true')
input_kwargs = vars(ap.parse_args())

machines_dct = cjson_load(input_kwargs['machines'])

# Trying to automatically determine the supercomputer if not set
if not input_kwargs['machine']:
    hostname = socket.gethostname()
    input_kwargs['machine'] = 'local'
    for mach_k, mach_v in machines_dct.items():
        if mach_v['hostname'] in hostname:
            input_kwargs['machine'] = mach_k
    print('Machine is automatically determined as: ' + input_kwargs['machine'])
    print('In case of error, specify the machine explicitly using -c key and edit the machines.cjson accordingly')

# Assigning input arguments to variables, minding the defaults 
# recalc_fold = input_kwargs['recalc_folder']
# poscars_fname = input_kwargs['poscars_file']
db_file = input_kwargs['database_file']
# maxcalcs = int(input_kwargs['maxcalcs'])
# outfilepattern = input_kwargs['outfile_pattern']


while not os.path.isfile('DONE') and not os.path.isfile('STOP'):
    # Checking if database file exists:
    print('*' * 10 + '{:%Y-%m-%d %H:%M}'.format(datetime.now()) + '*' * 10 + '\n')
    if os.path.isfile(db_file):
        with open(db_file, 'rb') as db_fid:
            database = pickle.load(db_fid)
        database.update(scenarios=input_kwargs['scenarios'])
    else:
        database = GauCalcDB(machines_json=machines_dct,
                             **{k: input_kwargs[k] for k in ('scenarios', 'poscars_file', 'recalc_folder', 'maxcalcs',
                                                             'machine', 'outfile_pattern', 'min_mult')})
        database.update(scenarios=input_kwargs['scenarios'])
    database.submit_jobs(n_par_calcs=int(input_kwargs['maxcalcs']))
    database.get_stats(verb=True)
    database.dump(filename_pkl=input_kwargs['database_file'], filename_en='energies.txt')

    # sys.stdout.flush()
    sleep_sec = int(float(input_kwargs['sleep']) * 60)
    time.sleep(sleep_sec)

# sys.stdout.close()