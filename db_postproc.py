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
from src.ext_software_io import parse_gout
from ase.io import read, write
from ase import Atoms
from glob import glob


def process_db_folder(db_fold, element_nos, res_folder=None):

    """
    Utility analysing the folder containing database pkl files and returning
    1. dictionary {(composition,): {'old_ind': [old_indices],
                                 'energies': [energies],
                                 'ccdata': [cclib_datas],
                                 'taskname': [task_names]}}  # for retreiving old indices

    2. list of information of lowest energy structures for each composition.
    Each element of the list is of the following format:
    [[composition], energy, gap, is_switched]
    Example for C6H12, which has lowest structure different from the one obtained from :


    @param db_fold: string, a folder where pkl's with GauCalcDB's are stored
    @param element_nos: tuple with periodic numbers of elements in the desired order, e.g. (6, 1) for C H
    @param res_folder: string, a folder where all desired data will be stored
    @return: :rtype: (dict, list)
    """

    if res_folder is None:
        resfolder = '.'
    if not os.path.exists(resfolder):
        os.makedirs(resfolder)

    list_fmt_best = []
    db_pkl_fnames = [str(x) for x in glob(db_fold + '/*.pkl')]
    for db_file in db_pkl_fnames:
        master_dict = {}
        with open(db_file, 'rb') as db_fid:
            res_poscars = db_file.split('.')[0] + '_res_POSCARS'
            database = pickle.load(db_fid)
            for task in database:
                try:
                    comp = tuple(np.count_nonzero(task.ccdata.atomnos == at_no) for at_no in element_nos)
                except AttributeError:
                    print(task.name + ' in ' + db_file + ' : job failed')
                    continue
                coords = task.ccdata.atomcoords[-1]
                coords -= np.sum(coords, axis=0)
                numbers = task.ccdata.atomnos
                cell = np.diag([np.max(coords[:, 0]) - np.max(coords[:, 0]) + 15,
                                np.max(coords[:, 1]) - np.max(coords[:, 1]) + 15,
                                np.max(coords[:, 2]) - np.max(coords[:, 2]) + 15])
                pbc = np.ones(3, dtype=bool)
                write(res_poscars, Atoms(positions=coords, numbers=numbers, pbc=pbc, cell=cell), append=True,
                      vasp5=True)
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
            homo_indices = val['ccdata'][lowest].homos
            if len(homo_indices) == 2:
                gap = np.min([val['ccdata'][lowest].moenergies[0][homo_indices[0] + 1],
                              val['ccdata'][lowest].moenergies[1][homo_indices[1] + 1]]) - \
                      np.max([val['ccdata'][lowest].moenergies[0][homo_indices[0]],
                              val['ccdata'][lowest].moenergies[1][homo_indices[1]]])
            else:
                gap = val['ccdata'][lowest].moenergies[0][homo_indices[0] + 1] - \
                      val['ccdata'][lowest].moenergies[0][homo_indices[0]]

            list_fmt_best.append([list(comp)] + [lowest_en] + [gap] + [changed])
            val['ccdata'][lowest].metadata['comments'] = 'E_tot = {:6.5f}'.format(lowest_en)
            val['ccdata'][lowest].writexyz(resfolder + '/' + val['taskname'][lowest] + '.xyz')

    # list_fmt_data = np.array(list_fmt_data)
    np.savetxt(resfolder + '/stats_np.txt', np.array(list_fmt_best)[:, :-1])
    list_fmt2table(np.array(list_fmt_best)[:, [0, 1, 2]], outfile=resfolder + '/en_table.txt')
    list_fmt2table(np.array(list_fmt_best)[:, [0, 1, 3]], outfile=resfolder + '/gap_table.txt')
    with open(resfolder + '/n_m_Enm_gap.txt', 'w') as stats_fid, open(resfolder + '/CH_gaps.txt', 'w') as gaps_fid:
        for dt in list_fmt_best:
            stats_fid.write('%d %d %.6f %.6f %s \n' % tuple(dt))
            gaps_fid.write('%2d\t%2d\t%.6f\n' % tuple(np.array(dt)[[0, 1, 3]]))

    return master_dict, list_fmt_best
