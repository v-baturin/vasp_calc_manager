#!python3
import sys
import argparse
import socket
import os
import cclib
from datetime import datetime
import pickle
import numpy as np
from tools_stability.aux_routines import list_fmt2table, table2list_fmt
import time
from src.common_tools import cjson_load
from src.tasks_database import GauCalcDB

element_nos = (6, 1)
resfolder = '/home/vsbat/SYNC/00__WORK/20191022_CH_project/res'
if not os.path.exists(resfolder):
    os.makedirs(resfolder)

list_fmt_data = []
db_pkl_fnames = ['ark_vb.pkl', 'ark_sl.pkl', 'ark_an.pkl', 'oleg_vb.pkl', 'oleg_sl.pkl', 'oleg_N.pkl', 'oleg_vb_pureC.pkl']
for db_file in db_pkl_fnames:
    master_dict = {}
    with open(db_file, 'rb') as db_fid:
        database = pickle.load(db_fid)
        for task in database:
            try:
                comp = tuple(np.count_nonzero(task.ccdata.atomnos == at_no) for at_no in element_nos)
            except AttributeError:
                print(task.name + ' in ' + db_file + ' : job failed')
                continue
            if comp in master_dict:
                master_dict[comp]['old_ind'].append(int(task.name.split('_')[-1]))
                master_dict[comp]['energies'].append(task.ccdata.scfenergies[-1])
                master_dict[comp]['ccdata'].append(task.ccdata)
                master_dict[comp]['taskname'].append(task.name)
            else:
                master_dict[comp] = {'old_ind': [int(task.name.split('_')[-1])],
                                     'energies': [task.ccdata.scfenergies[-1]],
                                     'ccdata': [task.ccdata],
                                     'taskname': [task.name]}
            print(task.name)

    # Processing and sorting
    for comp, val in master_dict.items():
        new_ind = np.argsort(val['energies'])
        lowest = new_ind[0]
        lowest_en = val['energies'][lowest]
        changed = (lowest != 0)
        list_fmt_data.append(list(comp) + [lowest_en] + [changed])
        val['ccdata'][lowest].metadata['comments'] = 'E_tot = {:6.5f}'.format(lowest_en)
        val['ccdata'][lowest].writexyz(resfolder + '/' + val['taskname'][lowest] + '.xyz')

# list_fmt_data = np.array(list_fmt_data)
np.savetxt(resfolder + '/stats_np.txt', np.array(list_fmt_data)[:,:-1])
list_fmt2table(np.array(list_fmt_data)[:,:-1], outfile=resfolder + '/en_table.txt')
with open(resfolder + '/stats.txt', 'w') as stats_fid:
    for dt in list_fmt_data:
        stats_fid.write('%d %d %.6f %s \n' % tuple(dt))


