#!python3
import os
import pickle
import numpy as np
from tools_stability.aux_routines import list_fmt2table, table2list_fmt
from ase.io import read, write
from ase import Atoms
from pathlib import Path
from cclib.parser.utils import PeriodicTable as pt
from genutils.abInitio_io_parse.gau_parse import get_gap, get_formula, get_kohn_sham
from matplotlib import pyplot as plt
from graph_utils.my_graphs import draw_spectrum
from graph_utils.formatting import cm2inch, set_ax_position_cm

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


def process_db_files(db_pkl_fnames,
                     element_nos=None, element_symbols=None,
                     res_folder=None,
                     n_isoms=1,
                     write_poscars = False,
                     write_xyz=False,
                     write_text = True,
                     first_connected_map=None):

    if element_nos and not element_symbols:
        element_symbols = tuple(element_labels[i] for i in element_nos)

    master_dict = db_files_to_dict(db_pkl_fnames,
                                   element_nos=element_nos,
                                   element_symbols=element_symbols)
    if write_poscars:
        write_db_to_poscars(master_dict,
                      res_folder + '/' + db_pkl_fnames[0].split('/')[-1].split('.')[0] + 'POSCARS')
    if write_xyz:
        write_xyz_of_n_lowest(master_dict, n_lowest=n_isoms, outdir=res_folder + '/xyz_files',
                              first_connected_map=first_connected_map)
    if write_text:
        best_list = list_of_bests(master_dict, first_connected_map=first_connected_map, lowest_inds_dict_out=False)
        write_txt_data(best_list, element_symbols=element_symbols, res_folder=res_folder)

    return master_dict




def db_files_to_dict(db_pkl_fnames,
                     element_nos=None, element_symbols=None,
                     res_folder=None,
                     skipnangap=False):

    """
    Utility analysing the folder containing database pkl files (recursively) and returning
    1. dictionary {(composition,): {'old_ind': [old_indices],
                                 'energies': [energies],
                                 'ccdata': [cclib_datas],
                                 'taskname': [task_names]}}  # for retreiving old indices



    @param db_pkl_fnames: list of paths to GauCalcDB's, or one path as a string
    @param element_nos: tuple with periodic numbers of elements in the desired order, e.g. (6, 1) for C H

    """
    if isinstance(db_pkl_fnames, str):
        db_pkl_fnames = [db_pkl_fnames]

    elif element_symbols and not element_nos:
        element_nos = tuple(ptable.number[sym] for sym in element_symbols)
    if res_folder is None:
        res_folder = '.'
    if not os.path.exists(res_folder):
        os.makedirs(res_folder)

    master_dict = {}
    for db_file in db_pkl_fnames:
        with open(db_file, 'rb') as db_fid:
            database = pickle.load(db_fid)
            for task in database:
                try:
                    comp = tuple(np.count_nonzero(task.ccdata.atomnos == at_no) for at_no in element_nos)
                except AttributeError:
                    print(task.name + ' in ' + db_file + ' : job failed')
                    continue
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
    return master_dict

def write_db_to_poscars(master_dict, poscar_fname, vacuumsize=15., append=True):
    sorted_cmp = get_sorted_compositions(master_dict)
    for cmp in sorted_cmp:
        comp = tuple(cmp)
        for ccd_i in master_dict[comp]['ccdata']:
            coords = ccd_i.atomcoords[-1]
            coords -= (np.sum(coords, axis=0)) / len(coords)
            numbers = ccd_i.atomnos
            cell = np.diag(np.max(coords, axis=0) - np.min(coords, axis=0) + vacuumsize)
            pbc = np.ones(3, dtype=bool)
            write(poscar_fname, Atoms(positions=coords, numbers=numbers, pbc=pbc, cell=cell), append=append,
                  vasp5=True)


def write_xyz_of_n_lowest(master_dict, n_lowest, outdir='xyz_files', first_connected_map=None):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    best_list, best_idcs = list_of_bests(master_dict, first_connected_map=first_connected_map, n_lowest=n_lowest, lowest_inds_dict_out=True)
    for cmp, val in master_dict.items():
        for i, lowest_k in enumerate(best_idcs[cmp]):
            val['ccdata'][lowest_k].metadata['comments'] = [get_formula(val['ccdata'][lowest_k]) +
                                                            ' dE = {:6.5f}'.format(val['energies'][lowest_k] -
                                                                                  val['energies'][best_idcs[cmp][0]])]
            val['ccdata'][lowest_k].writexyz(outdir + '/' + val['taskname'][lowest_k] + '_g' +
                                             str(i) + '.xyz')


def list_of_bests(master_dict, first_connected_map=None, n_lowest=1, lowest_inds_dict_out=True):
    """
    @param first_connected_map: BINARY SYSTEMS ONLY! (SiH, PdBi etc.)  2d array of numbers of lowest fully-connected
                                 isomers
    @return: :rtype: (dict, list)
    Return
    list of information of lowest energy structures for each composition.
    Each element of the list is of the following format:
    [[composition], energy, gap, is_switched]
    Example for C6H12, which has lowest structure different from the one obtained from :
    [[6, 1], -2255.445, 3.5, True]
    """

    list_fmt_best = []
    n_lowest_dict = {}

    sorted_cmp = get_sorted_compositions(master_dict)

    for cmp in sorted_cmp:
        comp = tuple(cmp)
        val = master_dict[comp]
        if first_connected_map is not None:
            lowest_connected = int(first_connected_map[comp])
        else:
            lowest_connected = 0
        new_ind = np.argsort(val['energies'])
        n_lowest_inds = new_ind[np.arange(lowest_connected, min(lowest_connected + n_lowest, len(new_ind)))]
        assert len(n_lowest_inds) > 0, "Invalid first_connected_map"
        lowest_en = val['energies'][n_lowest_inds[0]]
        changed = (n_lowest_inds[0] != 0)
        gap = val['gap'][n_lowest_inds[0]]
        list_fmt_best.append([list(comp)] + [lowest_en] + [gap] + [changed] + [n_lowest_inds[0]])
        if lowest_inds_dict_out:
            n_lowest_dict[comp] = n_lowest_inds
    if lowest_inds_dict_out:
        return list_fmt_best, n_lowest_dict
    else:
        return list_fmt_best


def write_txt_data(list_fmt_best, element_symbols, res_folder):
    n_el = len(list_fmt_best[0][0])
    if n_el == 2:
        flatten_list = [x[0] + x[1:-1] for x in list_fmt_best]
        list_fmt2table(np.array(flatten_list)[:, [0, 1, 2]], outfile=res_folder + '/en_table.txt')
        list_fmt2table(np.array(flatten_list)[:, [0, 1, 3]], outfile=res_folder + '/gap_table.txt')
    with open(res_folder + '/n_m_Enm_gap_2.txt', 'w') as stats_fid, open(res_folder + '/all_gaps.txt', 'w') as gaps_fid:
        stats_fid.write(('%s\t' * n_el) % element_symbols + 'Etot, eV\tGap, ev\tinitial ind\n')
        gaps_fid.write(('%s\t' * n_el) % element_symbols + '\tGap, ev\n')
        for dt in list_fmt_best:
            stats_fid.write(('%d\t' * n_el + '%.6f\t%.6f\t%s\t%d\n') % tuple(dt[0] + dt[1:]))
            gaps_fid.write(('%d\t' * n_el + '%.6f\n') % (tuple(dt[0]) + (dt[2],)))


def get_sorted_compositions(master_dict):
    """
    Utility for getting sorted list of compositions, where latest index is sorted first
    Example: dictionary consists of compositions [1,3], [2,4], [1,2], [2,1]
    Output: [[1,2], [1,3], [2,1], [2,4]]
    @param master_dict: master-dictionary extracted from database(s)
    @return: list of lists
    """
    composition_array = np.array([list(key) for key in master_dict.keys()])
    composition_array = composition_array[composition_array[:, -1].argsort()]
    for ind in range(-2, -composition_array.shape[1] - 1, -1):
        composition_array = composition_array[composition_array[:, ind].argsort(kind='mergesort')]
    return composition_array


def plot_el_spectra_binary(master_dict, element_symbols, savefiles=True, save_folder='spectra_graphs', groupped=True,
                           format='png', norm_by_natom=False, **drawspectrum_kwargs):

    # plt.rcParams['xtick.major.pad'] = '0'
    # plt.rcParams['ytick.major.pad'] = '0.'
    plt.rcParams['xtick.labelsize'] = 7
    plt.rcParams['ytick.labelsize'] = 7

    if savefiles and not os.path.isdir(save_folder):
        os.makedirs(save_folder)

    # master_dict =  db_files_to_dict(db_pkl_fnames, element_symbols=element_symbols)
    composition_array = get_sorted_compositions(master_dict)

    figno = 0
    for first_el, cnt in np.column_stack(np.unique(composition_array[:, 0], return_counts=True)):
        if groupped:
            figno += 1
            plt.figure()
            fig, axs = plt.subplots(cnt, 1, sharex=True)
            plt.setp(axs, yticks=np.array([-0.5, 0, 0.5]))
            fig.subplots_adjust(hspace=0)
            fig.set_size_inches(cm2inch(8, 1.5*cnt))
        for i, composition in enumerate(composition_array[composition_array[:,0] == first_el]):
            db_key = tuple(composition)
            db_val = master_dict[db_key]
            lowest_isom_idx = np.argsort(db_val['energies'])[0]
            cc_data = db_val['ccdata'][lowest_isom_idx]
            up_dos, dn_dos, fermi = get_kohn_sham(cc_data)

            up_dos -= fermi
            if cc_data.nelectrons % 2:
                dn_dos -= fermi
            else:
                dn_dos = None

            if groupped:
                base = axs[i]
            else:
                fig = plt.figure()
                fig.set_size_inches(cm2inch(8, 3))
                base = plt

            normalization = np.sum(composition) if norm_by_natom else True

            draw_spectrum(specup=up_dos, specdn=dn_dos, e_fermi=0, shareax=base,
                          label=get_formula(cc_data), sigma=0.05, normalization=normalization, **drawspectrum_kwargs)
            plt.xticks(np.arange(*drawspectrum_kwargs['span'], 1))
            if not groupped:
                ax = plt.gca()
                set_ax_position_cm(ax, [0.8, 0.8, 7, 2])
                ax.tick_params(bottom=True, top=True, direction='in')
                plt.savefig(save_folder + '/' + get_formula(cc_data) + '.' + format, dpi=600)
                plt.close()

        if groupped:
            plt.savefig(save_folder + '/' + element_symbols[0] + str(first_el) + '.' + format, dpi=600)
            plt.close()


    # compositions_dict = dict(zip(element_symbols, [[]] * len(element_symbols)))
    # for key, val in master_dict:
    #     fla = get_formula(val[0]['ccdata'], extended_out=True)['fdict']
    #     for el_sym in element_symbols:
    #         compositions_dict[el_sym].append(fla[el_sym])



if __name__ == '__main__':
    pass
