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

if len(sys.argv) == 1:
    sys.argv.extend(['-p', 'POSCARS_TEST',
                     '-f', 'RECALC',
                     '-d', 'database.pkl',
                     '-c', 'local', '--maxcalc', '3'])
print(sys.argv)  # debug

ap = argparse.ArgumentParser()
ap.add_argument("-c", "--computer", required=False, help='Which cluster')
ap.add_argument("-f", "--recalc_folder", required=False, help='Calc folder', default='RECALC')
ap.add_argument("-p", "--poscars_file", required=False, help="POSCARS-file", default='POSCARS')
ap.add_argument("-d", "--database_file", required=False, help='DB file', default='database.pkl')
ap.add_argument("-m", "--maxcalcs", required=False, help='Maximum parallel calcs', default='1')
ap.add_argument("-o", "--outfile_pattern", required=False, help='Maximum parallel calcs', default='log')
ap.add_argument("-s", "--sleep", required=False, help='Sleep time, in minutes', default='5')
ap.add_argument("--machines", required=False, help='Machines cjson file', default='machines.cjson')
ap.add_argument("--scenarios", required=False, help='Scenarios cjson file', default='scenarios.cjson')
args = vars(ap.parse_args())

machines_dct = cjson_load(args['machines'])

# Trying to automatically determine the supercomputer if not set
if args['computer']:
    machine = args['computer']
else:
    hostname = socket.gethostname()
    machine = 'local'
    for mach_k, mach_v in machines_dct.items():
        if mach_v['hostname'] in hostname:
            machine = mach_k
    print('Machine is automatically determined as: ' + machine)
    print('In case of error, specify the machine explicitly using -c key and edit the machines.cjson accordingly')

# Assigning input arguments to variables, minding the defaults 
recalc_fold = args['recalc_folder']
poscars_fname = args['poscars_file']
db_file = args['database_file']
maxcalcs = int(args['maxcalcs'])
outfilepattern = args['outfile_pattern']
sleep_sec = int(float(args['sleep']) * 60)

while not os.path.isfile('DONE') and not os.path.isfile('STOP'):
    # Checking if database file exists:
    print('*' * 10 + '{:%Y-%m-%d %H:%M}'.format(datetime.now()) + '*' * 10 + '\n')
    if os.path.isfile(db_file):
        with open(db_file, 'rb') as db_fid:
            database = pickle.load(db_fid)
        database.update(scenarios=args['scenarios'])
    else:
        database = GauCalcDB(scenarios=args['scenarios'], poscars_fname=poscars_fname, destination_fold=recalc_fold,
                             maxiter=maxcalcs, computer=machine, machines_json=machines_dct, outpattern=outfilepattern)

    database.submit_jobs(n_par_calcs=maxcalcs)
    database.get_stats(verb=True)
    database.dump(filename_pkl=db_file, filename_en='energies.txt')

    # sys.stdout.flush()

    time.sleep(sleep_sec)

# sys.stdout.close()