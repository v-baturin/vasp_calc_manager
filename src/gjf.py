import re
import warnings
import json
import math
import numbers
import os.path
import copy
from src.common_tools import recursively_map_to_vals, recursive_map_to_keys
from src.ext_software_io import gauformat2dict, dict2gauformat, xyz_2atoms, atoms_2_str


class Gjf(dict):
    """
    Assuming all numerical options for all keywords to be integer
    """

    memfactors = {'k': 0.1, 'm': 1., 'g': 1000., 't': 1e6}  # memory is internally represented in MB

    @staticmethod
    def gjfstring2dict(block):

        curr_gjf = dict()

        # Command section
        command_sect = re.findall('(%.*)(?:#)', block, re.DOTALL)
        curr_gjf['command'] = dict()
        if len(command_sect) == 1:
            command_sect = re.sub('[\n]*%', ' ', command_sect[0]).strip()
            curr_gjf['command'] = gauformat2dict(command_sect, case_sensitive=True)
            if 'mem' in curr_gjf['command'].keys():
                memory = curr_gjf['command']['mem'].casefold()
                prefix = re.search(r'[^ .0-9]+', memory).group(0)[0]
                numeric = re.search(r'[.0-9]+', memory).group(0)
                curr_gjf['command']['mem'] = int(math.ceil(float(numeric) * Gjf.memfactors[prefix]))
        elif len(command_sect) > 1:
            warnings.warn('Corrupted gjf file (command section)')
            return

        # Route section
        route_sect = re.search(r'(#.*?)(?=^$)', block, re.DOTALL | re.MULTILINE).group(0).replace('\n', ' ')
        if not route_sect:
            warnings.warn('Wrong file format (route section)')
            return
        curr_gjf['route'] = gauformat2dict(route_sect, gau_route=True)

        # Job name section
        jobname = re.findall('(?:#.*?\n\n)(.*?)(?=\n\n)', block, re.MULTILINE| re.DOTALL)
        if bool(jobname):
            curr_gjf['jobname'] = jobname[0]

        # Charge-multiplicity section
        chmult = re.findall('(?:\n\n)(-?\d+[\t ]+\d+[\t ]*)$', block, re.MULTILINE| re.DOTALL)
        if chmult:
            curr_gjf['charge_mult'] = chmult[0]

        # Molecular structure section
        molstruct = re.findall(r'((?:^[\t ]*[A-Z][a-z]?[\t ]+(?:[-0-9.]+.*){3}\n)+)', block, re.MULTILINE)
        if molstruct:
            curr_gjf['molstruct'] = xyz_2atoms(molstruct[0])
            if len(molstruct) == 2:
                print('Two geometries encountered')
        else:
            molstruct = None


        # Lasts section - custom basis, connectivity, wfn-file, etc
        last_section = re.findall('(?:-?\d+\s+\d+\s*.*?\n\n)(.*)(?=\n\n)', block, re.MULTILINE| re.DOTALL)
        if last_section:
            curr_gjf['last_section'] = last_section[0]
        return curr_gjf

    @staticmethod
    def apply_incremental_corrector(corrector: '10+1<20', value=None, key=None):
        (x0, pm, dx, dummy, thresh) = re.findall('(-?[0-9]+)([+-])(-?[0-9]+)([<>])(-?[0-9]+)', corrector)[0]
        [x0, dx, thresh] = [int(x) for x in (x0, dx, thresh)]
        if value is None:
            return x0
        elif isinstance(value, numbers.Number):
            if pm == '+':
                if value + dx <= thresh:
                    return value + dx
                else:
                    warnings.warn('Attempted to exceed maximum value. Keeping ' + str(key) + ' = ' + str(thresh))
                    return  thresh
            elif pm == '-':
                if value - dx >= thresh:
                    return value - dx
                else:
                    warnings.warn('Attempted to set value lower than minimum. Keeping ' + str(key) + ' = '
                                  + str(thresh))
                    return thresh
            else:
                raise ValueError("Wrong format of corrector. Should be 'n+d<m' or 'n-d>m'")

    @staticmethod
    def getdefaults(x):
        return Gjf.apply_incremental_corrector(x) if isinstance(x, str) else x

    @staticmethod
    def remove_minus(x):
        return x[1:] if x[0] == '-' else x

    def __init__(self, gjfdata, name='', meta=None):
        """
        Gjf class can be instantiated by:
        1. dictionary
        2. '<filename>.json' string
        3. a literal string representation of a dictionary: {'command':{'chk':'file.chk'...}}
        4. '<filename>.gjf' string
        @param gjfdata:
        @param name: object's name
        @param meta: dummy container
        """
        dict.__init__(self)
        self.name = name
        if meta:
            self.meta = meta
        if isinstance(gjfdata, dict):
            self.update(gjfdata)
        elif isinstance(gjfdata, str):
            if '.json' in gjfdata:
                with open(gjfdata) as f:
                    self.update(json.load(f))
            elif '\n' in gjfdata:
                self.update(Gjf.gjfstring2dict(gjfdata))
            else:
                if not os.path.isfile(gjfdata):
                    raise FileNotFoundError('File ' + gjfdata + " doesn't exist")
                print('Parsing ' + gjfdata)
                with open(gjfdata) as f:
                    block = ''.join(f)
                self.update(Gjf.gjfstring2dict(block))
                self.name = gjfdata.split('/')[-1]  # filename
                
    def copy(self):
        newone = Gjf(copy.deepcopy({k: v for k, v in self.items()}))
        newone.__dict__.update(self.__dict__)
        return newone

    def recurs_adjust(self, corrector, root_name=None):
        """
        Applies corrections specified in corrector to a Gjf-instance
        A corrector can be either dictionary (and Gjf-instance) or gjf-file
        Considering that 'key' is a name of gaussian KEYWORD and 'val' is a gaussian OPTION,
        The following behavior is specified:
            'key' in corrector.keys() AND 'key' in original.keys():
                val_corrector is dict and val_original is dict:
                    val_result = recurs_adjust(val_corrector, and val_original)
                (val_corrector is dict and val_original is either int or None) OR val_corrector is int:
                    val_result = val_corrector
                val_corrector is of form 'd[+/-]i[<>]m', d,i,m - int numbers (e.g. '100+50<1000' or '150-10>-500'):
                    val_original is either None or dict:
                        val_result = d
                    val_original is int:
                        val result = val_original[+/-]i if val_original[+/-]i [<>] m else: ERROR
            'key=' in corrector.keys() OR 'key' is NOT in original.keys()
                val_result = val_corrector (with all lower 'd[+/-]i[<>]m' turned to d)
            '-key' in corrector.keys()
                del original['key']
        @param corrector: correcting dictionary or filename
        @param root_name: name of the object (if None, .name of Gjf-object is taken if possible)
        @return: corrected self
        """
        if isinstance(corrector, str):
            corrector = Gjf(corrector)

        if root_name is None:
            if isinstance(self, Gjf):
                root_name = self.name
            elif isinstance(self, dict):
                root_name = ''

        for corr_key, corr_val in corrector.items():
            if corr_key in self.keys():
                orig_val = self[corr_key]
                if isinstance(corr_val, dict) and isinstance(orig_val, dict):
                    Gjf.recurs_adjust(orig_val, corr_val, root_name)
                elif isinstance(corr_val, str):
                    if re.match(r'(-?\d+\+-?\d+<-?\d+)|(-?\d+--?\d+>-?\d+)', corr_val):
                        self[corr_key] = Gjf.apply_incremental_corrector(corr_val, value=orig_val, key=corr_key)
                    else:
                        self[corr_key] = corr_val  # For non-route sections
                        warnings.warn('In ' + str(root_name) + ' ' + corr_key + ': '
                                      + str(orig_val) + ' is replaced by ' + str(corr_val))
                elif corr_val is None:
                    pass
                else:
                    step = recursive_map_to_keys(Gjf.remove_minus, corr_val)
                    self[corr_key] = recursively_map_to_vals(Gjf.getdefaults, step)
                    if not isinstance(orig_val, str):
                        warnings.warn('In ' + str(root_name) + ' ' + corr_key + ': '
                                      + str(orig_val) + ' is replaced by ' + str(corr_val))
            elif corr_key[0] == '-':
                if corr_key[1:] in self.keys():
                    del self[corr_key[1:]]
            else:
                if corr_key[-1] == '=':
                    new_corr_key = corr_key[:-1]
                else:
                    new_corr_key = corr_key

                if isinstance(corr_val, dict):
                    step = recursive_map_to_keys(Gjf.remove_minus, corr_val)
                    self[new_corr_key] = recursively_map_to_vals(Gjf.getdefaults, step)
                elif isinstance(corr_val, str) and re.match(r'(-?\d+\+-?\d+<-?\d+)|(-?\d+--?\d+>-?\d+)', corr_val):
                    self[new_corr_key] = Gjf.getdefaults(corr_val)
                else:
                    self[new_corr_key] = corr_val

    def save(self, fname, openmode='w'):

        with open(fname, mode=openmode) as fid:

            # Command section
            if 'command' in self.keys():
                for key, val in self['command'].items():
                    if key == 'mem':
                        fid.write('%' + key + '=%sMB\n' % val)
                    else:
                        fid.write('%' + key + '=%s\n' % val)

            # Route section
            route_str = dict2gauformat(self['route'], gau_route=True)
            fid.write(route_str + '\n')

            if route_str != '#p restart':
                # Job name
                if 'jobname' in self.keys():
                    fid.write('\n' + self['jobname'].strip() + '\n')

                # Charge multiplicity
                if 'charge_mult' in self.keys():
                    fid.write('\n' + self['charge_mult'].strip() + '\n')

                # Molecular structure
                if 'molstruct' in self.keys():
                    fid.write(atoms_2_str(self['molstruct']) + '\n')


                if 'last_section' in self.keys():
                    fid.write('\n' + self['last_section'] + '\n')

            fid.write('\n')


if __name__ == '__main__':
    # # fname = '../example.gjf'
    # # gj1 = Gjf(fname)
    # # print(gj1)
    # fname = '../testfolder/complicated.gjf'
    # mooname = '/home/vsbat/SYNC/00_Current_PyWork/refine_with_gaussian/testfolder/Mo10O14_orb.gjf'
    # gj2 = Gjf(mooname)
    # print(gj2)
    # test_corr = {'route': {'output': {'wfn': 5}}}
    # gj2.recurs_adjust(test_corr)
    # print(gj2)
    # gj2.save('nuka.gjf')
    # 
    # testorig = {'k1': {'kk': None, 'k2': 3, 'k3': {'k4': {'hello': None}, 'k5': 10}}}
    # # testcorr = {'k1': {'kk': 4,'k2': '4+2<6', 'k3=': {'k4': None, '-k5': None, 'k6': {'k7':None, 'k8': '5+2<15'}, 'k9': '7-1>3'}}}
    # testcorr = {'-k1': None}
    # gj3 = Gjf(testorig)
    # gj3.recurs_adjust(testcorr)
    # print(gj3)

    # gjf_fname = '../gjf_templates/CH.gjf'
    # gj4 = Gjf(gjf_fname)
    # print(gj4['command'])

    gjf_fname = '/home/vsbat/mnt/arkuda_VBATURIN/PROJECT_PdBi/TS/TS_attempt2/Kurzina/TS_for_1st_step/TS_for_1st_step.gjf'
    gj5 = Gjf(gjf_fname)

    corrector = {"route=": {"restart": None}}
    gj5.recurs_adjust(corrector)
    gj5.save('/home/vsbat/Desktop/testgj')
    print(gj5['command'])
