#!python3
import sys
import argparse
import socket
import os
from datetime import datetime
import pickle
import time
from src.common_tools import cjson_load, add_index
from src.tasks_database import GauCalcDB

start_time = time.time()

DEFAULT_RES_FOLDER = 'results'

ap = argparse.ArgumentParser()
ap.add_argument("-c", "--machine", required=False, help='Which cluster')
ap.add_argument("-f", "--recalc_folder", required=False, help='Calc folder')
ap.add_argument("-g", "--geoms_file", required=False, help="POSCARS-file or xyz-batch", default='POSCARS')
ap.add_argument("-m", "--maxcalcs", required=False, help='Maximum parallel calcs', default='1')
ap.add_argument("-o", "--outfile_pattern", required=False, help='Output file pattern', default='log')
ap.add_argument("-s", "--sleep", required=False, help='Sleep time, in minutes', default='5')
ap.add_argument("--maxiter", required=False, help='Maximum number of iteration', default='5')
ap.add_argument("--machines", required=False, help='Machines cjson file', default='machines.cjson')
ap.add_argument("--scenarios", required=False, help='Scenarios cjson file', default='scenarios.cjson')
ap.add_argument("--min_mult", required=False, help='Use minimal multiplicity', action='store_true')
ap.add_argument("--exec_time", required=False, help='Execution time, mins',  default='300')
input_kwargs = vars(ap.parse_args())

if input_kwargs['recalc_folder'] is None:
    input_kwargs['recalc_folder'] = add_index(DEFAULT_RES_FOLDER)

log_path = os.path.join(input_kwargs['recalc_folder'], 'log')
os.makedirs(input_kwargs['recalc_folder'], exist_ok=True)
sys.stdout = open(log_path, 'a')

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
db_file = os.path.join(input_kwargs['recalc_folder'], 'database.pkl')

if os.path.isfile(db_file):
    with open(db_file, 'rb') as db_fid:
        database = pickle.load(db_fid)
    db_file = db_file.split('/')[-1]
else:
    database = GauCalcDB(machines_json=machines_dct,
                         **{k: input_kwargs[k] for k in ('scenarios', 'geoms_file', 'recalc_folder', 'maxiter',
                                                         'machine', 'outfile_pattern', 'min_mult')})


while not os.path.isfile(os.path.join(input_kwargs['recalc_folder'], 'DONE')) and \
        not os.path.isfile(os.path.join(input_kwargs['recalc_folder'], 'STOP')) and \
        time.time() - start_time <= input_kwargs['exec_time']:
    database.update(scenarios=input_kwargs['scenarios'])
    print('*' * 10 + '{:%Y-%m-%d %H:%M}'.format(datetime.now()) + '*' * 10 + '\n')
    database.submit_jobs(n_par_calcs=int(input_kwargs['maxcalcs']))
    database.get_stats(verb=True)
    database.dump(filename_pkl=db_file, filename_en='energies.txt')
    sys.stdout.flush()
    sleep_sec = int(float(input_kwargs['sleep']) * 60)
    time.sleep(sleep_sec)

# sys.stdout.close()