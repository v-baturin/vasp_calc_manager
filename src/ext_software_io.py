import re
import json
import numpy as np
from copy import deepcopy
from ase import Atoms
from ase.io import read as ase_read
from cclib.io import ccread
from .common_tools import recursive_map_to_keys, recursively_map_to_vals, try_numerize_string


def read_poscars(poscars_fname):
    """
    Reads POSCARS multiple images file into a list of Atoms objects
    @param poscars_fname: POSCARS filename
    @return: list of Atoms objects
    """
    poscars_list = list()
    with open(poscars_fname) as f:
        while True:
            try:
                poscars_list.append(ase_read(f, format='vasp'))
            except IndexError:
                break
    return poscars_list


def gauformat2dict(string: str, gau_route=False, case_sensitive=False):
    """
    Converts a string from Gaussian input-file into a dict variable
    @param string: string of a specific card within Gaussian gjf file
    @param gau_route: boolean flag for Route section
    @param case_sensitive:
    @return: dictionary {k1:v1,k2:{k3:v3,k4:v4}...}
    """
    out_dict = dict()
    if not case_sensitive:
        string = string.casefold()
    string = string.strip()
    string = re.sub(' +', ' ', string)  # Removing redundant spaces
    string = re.sub('\s*=\s*', '=', string)  # Removing spaces around '='
    string = re.sub('\s*\)', ')', string)  # Removing spaces before ')'
    string = re.sub('\s*\(\s*', '(', string)  # Removing spaces before ')'
    string = re.sub('\s*,\s*', ',', string)  # Removing spaces around '('
    if gau_route:  # Special tratment for Route ('# ...') section
        string = re.sub(r'^#\w? +', '', string)  # removing standard beginning with trailing spaces '#p  '
        string = re.sub(r'^#', '', string)  # removing single hash before approach section like e.g. in #b3lyp/cam
        list_from_string = string.split()
        for k, item in enumerate(list_from_string):
            if re.match('[^/ ()]+/.*', item):
                approach_basis = item
                del list_from_string[k]
            else:
                approach_basis = ''
        string = ' '.join(list_from_string)
    string = re.sub(r'([^=])\(', r'\1=(', string)  # Inserting '=' sign for homogeneity. (par(vals) -> par=(vals))
    list_from_string = string.split()
    step = ','.join(list_from_string)
    step = re.sub(r'([^ \(\),=]+)', r"'\1'", step)  # Encloses all parameters in apostrophes
    step = re.sub('(^|[,\(])(\'[^ \(\),=]+\')(?=([,\)]|$))', r"\1\2:null",
                  step)  # insertion of colons for params without options
    step = re.sub('=\(', ':{', step)
    step = re.sub('\)', '}', step)
    finalstr = '{' + re.sub('=', ':', step) + '}'
    out_dict.update(json.loads(finalstr.replace('\'', '\"')))
    # Now we need to bring all keys to lower case and make numerical values really numerical:
    lower = lambda line: line.casefold()
    out_dict = recursive_map_to_keys(lower, recursively_map_to_vals(try_numerize_string, out_dict))
    # replace all str values s by single-item dictionaries {s: None}:
    if gau_route:
        todict = lambda v: {v: None} if isinstance(v, str) else v
        out_dict = recursively_map_to_vals(todict, out_dict)
        out_dict['approach_basis'] = approach_basis
    return out_dict


def dict2gauformat(dct_param: dict, gau_route=False):
    txt = ''
    dct = deepcopy(dct_param)
    if gau_route:
        if dct_param == {'restart': None}:
            return '#p restart'
        txt = txt + '#p ' + dct['approach_basis']
        del dct['approach_basis']
    for upper_key, upper_val in dct.items():
        txt += (' ' + upper_key)
        if isinstance(upper_val, dict):
            if len(upper_val) == 1 and next(iter(upper_val.values())) is None:
                txt += ('=' + next(iter(upper_val.keys())))
            elif upper_key == 'iop':
                txt += ('(' + recursive_print(upper_val) + ')')
            else:
                txt += ('=(' + recursive_print(upper_val) + ')')
        elif bool(upper_val):
            txt += ('=' + str(upper_val))
    return txt.replace(',)', ')')


def recursive_print(dct: dict):
    txt1 = ''
    for k, v in dct.items():
        if isinstance(v, dict):
            if len(v) == 1 and next(iter(v.values())) is None:
                txt1 += (k + '=' + next(iter(v.keys())))
            else:
                txt1 += (k + '=(' + recursive_print(v) + '),')
        elif bool(v):
            txt1 += (k + '=' + str(v) + ',')
        else:
            txt1 += (k + ',')
    return txt1.replace(',)', ')')


def xyz_2atoms(xyz):
    symbols = []
    positions = []
    closefile = False
    if isinstance(xyz, str):
        if '\n' in xyz:
            xyz = xyz.split('\n')
        else:
            closefile = True
            xyz = open(xyz)

    for item in xyz:
        if re.match(r'(^\s*[A-Z][a-z]?\s+(?:[-0-9.]+.*){3})', item):
            tokens = item.split()
            symbol = tokens[0]
            pos = list(map(float, tokens[1:4]))
            symbols.append(symbol)
            positions.append(pos)
    atoms = Atoms(symbols, positions)
    if closefile:
        xyz.close()
    return atoms


def gaussian_in_2_Atoms(fd):
    """
    Borrowed from ASE package
    @param fd:
    @return:
    """
    re_chgmult = re.compile(r'^\s*[+-]?\d+(?:,\s*|\s+)[+-]?\d+\s*$')

    symbols = []
    positions = []
    pbc = np.zeros(3, dtype=bool)
    cell = np.zeros((3, 3))
    npbc = 0
    # We're looking for charge and multiplicity
    for line in fd:
        if re_chgmult.match(line) is not None:
            tokens = fd.readline().split()
            while tokens:
                symbol = tokens[0]
                pos = list(map(float, tokens[1:4]))
                if symbol.upper() == 'TV':
                    pbc[npbc] = True
                    cell[npbc] = pos
                    npbc += 1
                else:
                    symbols.append(symbol)
                    positions.append(pos)
                tokens = fd.readline().split()
            atoms = Atoms(symbols, positions, pbc=pbc, cell=cell)
            return atoms


def atoms_2_str(atms):
    output = []
    for at in atms:
        output.append('{:<10s}{:20.10f}{:20.10f}{:20.10f}'.format(at.symbol, *at.position))
    return '\n'.join(output)


def parse_gout(logfile):
    print('parsing: ' + logfile)
    ccdata = ccread(logfile)
    try:
        last_positions = ccdata.atomcoords[-1]
    except:
        return None, None
    numbers = ccdata.atomnos
    return ccdata, Atoms(positions=last_positions, numbers=numbers)

# def get_options(string, option):
#     """
#     :param string: 'option=val option1=(val1, val2) option3(val3, val4)
#     :param option: one of the options
#     :return: list of values for given option
#     """
#     # Preconditioning the string:
#     string = re.sub(' +', ' ', string)
#     string = re.sub('\s*=\s*', '=', string)
#     string = re.sub('\s*\)', ')', string)
#     string = re.sub('\(\s*', '(', string)
#     string = re.sub('\s*,\s*', ',', string)
#
#     # Analyzing the string
#     if option.casefold() not in string.casefold():
#         return None
#     option_opts_findres = re.findall('(?:' + option + '[=\(]+)([^ \)\n]+)', string, re.IGNORECASE)
#
#     if len(option_opts_findres) == 0:
#         return []
#     else:
#         option_opts = ','.join(option_opts_findres).split(',')
#     return [j.strip() for j in option_opts]

if __name__ == '__main__':
    upper = lambda x: x.casefold()
    testdct = {'Mem': 10, 'HaHa': {'Lalal': 'AAA', 'uUuU': 777}}
    resdct = recursive_map_to_keys(upper, testdct)
    print(resdct)


    def try_numerize_string(string):
        if not isinstance(string, str):
            return string
        else:
            string = string.strip()
        if re.match('-?[0-9]+$', string):
            return int(string)
        else:
            try:
                return float(string)
            except:
                return string


    testdct2 = {'a': 'Hello', 'b': '150.00', 'c': 10, 'd': '10', 'e': '1e5'}
    resdct = recursively_map_to_vals(try_numerize_string, testdct2)
    print(resdct)

    g_in_fname = '/home/vsbat/SYNC/1_TALKS/20200604_Sk_lab_TALK/TXM-CAT_CHCl3_Step1.gjf'
    with open(g_in_fname) as f:
        ats = gaussian_in_2_Atoms(f)
    xyz_at = xyz_2atoms(g_in_fname)

    bad_g_out = '/home/vsbat/SYNC/00_Current_PyWork/refine_with_gaussian/testfolder/log737771.8835'
    # output_g = '/home/vsbat/SYNC/00_Current_PyWork/refine_with_gaussian/testfolder/log737771.7237'
    ccdata = ccread(bad_g_out)


    pass
    # poscar_fname = '/home/vsbat/SYNC/00_Current_PyWork/20200319_Catalysis_project/CuAu/CuAu_varcomp_6_10_6_10/results2/POSCAR'
    # poscar_atoms = ase_read(poscar_fname)
