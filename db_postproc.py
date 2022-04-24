#!python3
import os
import pickle
import numpy as np
from tools_stability.aux_routines import list_fmt2table, table2list_fmt
from ase.io import read, write
from ase import Atoms
# from glob import glob
from pathlib import Path
from cclib.parser.utils import PeriodicTable as pt
from genutils.abInitio_io_parse.gau_parse import get_gap

ptable = pt()
element_labels = np.array(ptable.element[:])

def process_db_folder(db_fold, *args, **kwargs):
    """
    @param db_fold: string, a folder where pkl's with GauCalcDB's are stored
    @param args: arguments to be passed to process_db_files
    @param kwargs: keyword arguments to be passed to process_db_files
    @return:
    """
    db_pkl_fnames = [str(x) for x in Path(db_fold).rglob('*.pkl')]
    return process_db_files(db_pkl_fnames, *args, **kwargs)


def process_db_files(db_pkl_fnames, element_nos, res_folder=None, n_isoms=1, write_xyz=False, skipnangap=False, first_connected_map=None):

    """
    Utility analysing the folder containing database pkl files (recursively) and returning
    1. dictionary {(composition,): {'old_ind': [old_indices],
                                 'energies': [energies],
                                 'ccdata': [cclib_datas],
                                 'taskname': [task_names]}}  # for retreiving old indices

    2. list of information of lowest energy structures for each composition.
    Each element of the list is of the following format:
    [[composition], energy, gap, is_switched]
    Example for C6H12, which has lowest structure different from the one obtained from :
    [[6, 1], -2255.445, 3.5, True]

    @param db_pkl_fnames: list of paths to GauCalcDB's, or one path as a string
    @param element_nos: tuple with periodic numbers of elements in the desired order, e.g. (6, 1) for C H
    @param res_folder: string, a folder where all desired data will be stored
    @param first_connected_map: BINARY SYSTEMS ONLY! (SiH, PdBi etc.)  2d array of numbers of lowest fully-connected
                                 isomers
    @return: :rtype: (dict, list)
    """
    VACUUM_SIZE = 15
    if isinstance(db_pkl_fnames, str):
        db_pkl_fnames = [db_pkl_fnames]

    element_symbols = tuple(element_labels[i] for i in element_nos)
    if res_folder is None:
        res_folder = '.'
    if not os.path.exists(res_folder):
        os.makedirs(res_folder)

    list_fmt_best = []
    master_dict = {}
    for db_file in db_pkl_fnames:
        with open(db_file, 'rb') as db_fid:
            db_fname = db_file.split('/')[-1].split('.')[0]
            res_poscars = db_fname + '_res_POSCARS'
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
                cell = np.diag(np.max(coords, axis=1) - np.min(coords, axis=1) + VACUUM_SIZE)
                pbc = np.ones(3, dtype=bool)
                write(res_poscars, Atoms(positions=coords, numbers=numbers, pbc=pbc, cell=cell), append=True,
                      vasp5=True)
                gap = get_gap(task.ccdata)
                if (skipnangap or np.isfinite(gap)):
                    if comp in master_dict:
                        master_dict[comp]['old_ind'].append(int(task.name.split('_')[-1]))
                        master_dict[comp]['energies'].append(task.ccdata.scfenergies[-1])
                        master_dict[comp]['ccdata'].append(task.ccdata)
                        master_dict[comp]['taskname'].append(task.name)
                        master_dict[comp]['where'].append(db_file.split('/')[-1])
                        master_dict[comp]['gap'].append(gap)
                    else:
                        master_dict[comp] = {'old_ind': [int(task.name.split('_')[-1])],
                                             'energies': [task.ccdata.scfenergies[-1]],
                                             'ccdata': [task.ccdata],
                                             'taskname': [task.name],
                                             'where': [db_file.split('/')[-1]],
                                             'gap':[gap]}
                # print(task.name)

    # Processing and sorting
    for comp, val in master_dict.items():
        if first_connected_map is not None:
            lowest_connected = int(first_connected_map[comp])
        else:
            lowest_connected = 0
        print(comp, len(val['energies']))
        new_ind = np.argsort(val['energies'])
        n_lowest_inds = new_ind[np.arange(lowest_connected, min(lowest_connected + n_isoms, len(new_ind)))]
        assert len(n_lowest_inds) > 0, "Invalid first_connected_map"
        lowest_en = val['energies'][n_lowest_inds[0]]
        changed = (n_lowest_inds[0] != 0)
        gap = val['gap'][n_lowest_inds[0]]
        list_fmt_best.append([list(comp)] + [lowest_en] + [gap] + [changed])
        if write_xyz:
            for i, lowest_k in enumerate(n_lowest_inds):
                val['ccdata'][lowest_k].metadata['comments'] = ['dE = {:6.5f}'.format(val['energies'][lowest_k] -
                                                                                      val['energies'][n_lowest_inds[0]])]
                val['ccdata'][lowest_k].writexyz(res_folder + '/' + val['taskname'][lowest_k] + '_g' +
                                                 str(i + lowest_connected) + '.xyz')

    # list_fmt_data = np.array(list_fmt_data)
    # with open(res_folder + '/stats_np.txt', 'w') as stats_
    #     np.savetxt(res_folder + '/stats_np.txt', np.array(list_fmt_best)[:, :-1])
    if len(element_nos) == 2:
        flatten_list = [x[0] + x[1:-1] for x in list_fmt_best]
        list_fmt2table(np.array(flatten_list)[:, [0, 1, 2]], outfile=res_folder + '/en_table.txt')
        list_fmt2table(np.array(flatten_list)[:, [0, 1, 3]], outfile=res_folder + '/gap_table.txt')
    with open(res_folder + '/n_m_Enm_gap.txt', 'w') as stats_fid, open(res_folder + '/all_gaps.txt', 'w') as gaps_fid:
        stats_fid.write(('%s\t' * len(element_nos)) % element_symbols + 'Etot, eV\tGap, ev\n')
        for dt in list_fmt_best:
            stats_fid.write(('%d\t'*len(element_nos) + '%.6f\t%.6f\t%s\n') % tuple(dt[0] + dt[1:]))
            gaps_fid.write(('%d\t'*len(element_nos) + '%.6f\n') % (tuple(dt[0]) + (dt[2],)))

    return master_dict, list_fmt_best


if __name__ == '__main__':
    pass
