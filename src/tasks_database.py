import re
import pickle
import os
import glob
from datetime import datetime
from os.path import isfile
from pathlib import Path
from .ext_software_io import read_poscars, parse_gout, dict2gauformat
from .gjf import Gjf
from .common_tools import add_index, sh_execute, checkflag, cjson_load, mk_new_dir

# termination_dict = {'init_gjf':  'initial_config.gjf',
#                     'normal':   {'flags' : 'Normal termination', 'tail': 1}, # To be modified for orb_coeff etc
#                     'restart':  {'corrector': {'route': {'scf': {'yqc': None}}}},
#                     'errors':   {'low_phys_mem': {'flags': 'could not allocate memory', 'tail': 20,
#                                                   'msg': 'Could not allocate memory',
#                                                   'corrector': {'command':{'mem': '4000-100<1000'}}},
#                                  'low_alloc_mem': {'flags': 'Out-of-memory error', 'tail': False,
#                                                    'msg': 'Insufficient memory allocated',
#                                                    'corrector': {'command':{'mem': '1000+1000<10000'}}}},
#                     'fails':    {'bad_input': {'flags': 'Error termination via.*l101', 'tail': 50,
#                                                'msg': 'Bad input file',
#                                                'corrector': {'route': {'scf': {'yqc': None}}}}}}


class GauTask:
    # machines = {'rurik': {'jobtemplate': 'job_templates/rurik.sh',
    #                       'runcmd': r'sbatch ', 'jobid': '(?<=job )\d+',
    #                       'queue': "squeue -u `whoami` | grep -Po '(^\s*)\d+' |grep -Po '\d+'"},
    #             'oleg': {'jobtemplate': 'job_templates/oleg.sh',
    #                      'runcmd': r'sbatch ', 'jobid': '(?<=job )\d+',
    #                      'queue': "squeue -u `whoami` | grep -Po '(^\s*)\d+' |grep -Po '\d+'"},
    #             'arkuda': {'jobtemplate': 'job_templates/arkuda.sh',
    #                        'runcmd': r'qsub ', 'jobid': '\d+(?=\.)', 'queue': "qstat | grep -Po '^\s*\d+'"},
    #             'local': {'jobtemplate': 'job_templates/local.sh',
    #                       'runcmd': r'', 'jobid': '\d+', 'queue': "ps -x | grep -Po '^\s*\d+'"}}

    def __init__(self, gjf=None, folder=None, status='P', k_iter=0, ccdata=None, job_id=None, name=None,
                 jobfile=None, machine='local', out_fname='log', machine_dict=None, chk_file=None):

        if not machine_dict:
            machine_dict = {'local': {'jobtemplate': 'job_templates/local.sh',
                                      'runcmd': r'', 'jobid': '\d+', 'queue': "ps -x | grep -Po '^\s*\d+'"}}

        self.job_id = job_id
        self.ccdata = ccdata
        self.k_iter = k_iter
        self.status = status
        self.gjf = gjf.copy()
        self.previousgjf = gjf.copy()
        if not chk_file:
            chk_file = name + '.chk'
        self.gjf.recurs_adjust({'command': {'chk=': chk_file}, 'jobname=': name})
        if 'molstruct' in gjf:
            self.old_coords = gjf['molstruct'].copy()
        self.folder = folder  # Full path
        self.name = name  # E.g. Cd2Se10_3
        self.jobfile = jobfile
        self.machine = machine
        self.machine_dict = machine_dict
        self.out = out_fname

    def __repr__(self):
        return '< Gaussian task ' + self.name + ' at ' + self.folder + '>'
    
    @property
    def gjffilename(self):
        return self.name + '.gjf'

    @property
    def full_gjf_path(self):
        return self.folder + '/' + self.gjffilename

    def isrunning(self):
        if self.job_id:
            cur_queue = sh_execute(self.machine_dict['queue'])
            return self.job_id in cur_queue
        else:
            return False

    def copyfiles(self):
        if not os.path.isdir(self.folder):
            os.makedirs(self.folder)
        self.gjf.save(self.full_gjf_path)
        if 'command' not in self.gjf.keys() or 'nprocshared' not in self.gjf['command']:
            self.gjf['command'] = {'nprocshared': ''}
        sh_execute("sed 's/@INPUT/" + self.gjffilename + "/' " + self.machine_dict['jobtemplate'] + " | " +
                   "sed 's/@OUTFILE/" + self.out + "/' | " +
                   "sed 's/@JOBNAME/" + self.name + "/' | " +
                   "sed 's/@NPROCSHARED/" + str(self.gjf['command']['nprocshared']) + "/' " +
                   " > " + self.folder + '/' + self.jobfile)
        sh_execute("chmod u+x " + self.folder + '/' + self.jobfile)

    def submit_job(self):
        print('Submitting ' + self.name)
        submit_cmd = 'cd ' + self.folder + ' && ' + self.machine_dict['runcmd'] + ' ./' + self.jobfile
        submit_out = sh_execute(submit_cmd, add_errors=True)
        print(submit_out)
        job_id = re.findall(self.machine_dict['jobid'], submit_out)[0]
        self.job_id = job_id
        print(job_id)


class GauCalcDB(list):

    def __init__(self, scenarios='scenarios.cjson', poscars_file=None, recalc_folder=None, maxiter=5,
                 machine='local', machines_json='machines.json', outfile_pattern='log', inpattern='*.gjf', min_mult=False):

        list.__init__(self)

        self.maxiter = int(maxiter)
        self.master_folder = recalc_folder

        scen_dct = cjson_load(scenarios)
        init_gjf = Gjf(scen_dct['init_gjf'])

        mach_dcts = cjson_load(machines_json)
        machine_dict = mach_dcts[machine]

        if poscars_file and isfile(poscars_file):

            poscars_list = read_poscars(poscars_file)

            for poscar in poscars_list:
                # Creating task entry
                formula = poscar.get_chemical_formula()
                init_gjf['molstruct'] = poscar
                if min_mult:
                    init_gjf['charge_mult'] = init_gjf['charge_mult'].split()[0] +\
                                              ' ' + str(sum(poscar.get_atomic_numbers()) % 2 + 1)
                curr_dest_folder = mk_new_dir(recalc_folder + '/' + formula + '/' + formula, zerobased=True)
                curr_label = curr_dest_folder.split('/')[-1]
                task = GauTask(gjf=init_gjf, folder=curr_dest_folder, name=curr_label, jobfile='jobfile.sh',
                               machine=machine, machine_dict=machine_dict,
                               out_fname=re.sub('.*\.', curr_label + '.', outfile_pattern))
                # task.copyfiles()  # Copying necessary files
                self.append(task)

        elif inpattern and outfile_pattern:

            all_logs = Path(recalc_folder).rglob(inpattern)
            n_logs = 0
            for k_log_file in all_logs:
                curr_folder_full = str(k_log_file.parent)
                curr_folder = str(k_log_file.parts[-2])
                work_gjf_name = sh_execute('ls ' + curr_folder_full + '/*.gjf').strip()

                jobfnamesearch = glob.glob(curr_folder_full + '/*.sh')
                if jobfnamesearch:
                    jobfname = jobfnamesearch[0].split('/')[-1]
                else:
                    jobfname = 'jobfile.sh'

                out_fname_search = glob.glob(curr_folder_full + '/' + outfile_pattern)
                if out_fname_search:
                    out_fname = out_fname_search[0].split('/')[-1]
                    newstatus = 'L'
                else:
                    out_fname = outfile_pattern.replace('*', 'log')
                    newstatus = 'P'

                chksearch = glob.glob(curr_folder_full + '/*.chk')
                chk_fname = chksearch[0].split('/')[-1] if chksearch else None

                gjf = Gjf(work_gjf_name)
                name = work_gjf_name.split('/')[-1].split('.')[0]
                task = GauTask(gjf=gjf, folder=curr_folder_full, name=name, jobfile=jobfname,
                               machine=machine, machine_dict=machine_dict, out_fname=out_fname, chk_file=chk_fname,
                               status=newstatus)
                self.append(task)
                n_logs += 1
            if n_logs == 0:
                raise Exception('Destination folder has no Gaussian calculations with in- and out-files of provided format!')
        else:
            raise Exception('Error: destination_folder and (poscars_file or outfile_pattern) must be provided! ')

    def __str__(self):
        return '<GauCalcDB: ' + str(len(self)) + ' GauTask objects in ' + self.master_folder + '>'

    def submit_jobs(self, n_par_calcs=3):
        submitted = new = 0
        for i, task in enumerate(self):
            if task.status == 'R':
                submitted += 1
            if task.status == 'P':
                task.copyfiles()
                try:
                    task.submit_job()
                except IndexError:
                    print('Failed to submit: ' + task.name)
                else:
                    self[i].status = 'R'
                    submitted += 1
                    new += 1
            if submitted >= n_par_calcs:
                break
        print('Submitted = ' + str(new) + '; Total running = ' + str(submitted))
        return new, submitted

    def update(self, scenarios='scenarios.cjson'):

        scenarios_dct = cjson_load(scenarios)

        for task in self:

            if not os.path.isfile(task.full_gjf_path):
                task.copyfiles()
                task.status = 'P'

            logfile = task.folder + '/' + task.out

            if task.k_iter >= self.maxiter:
                print(task.name + ': MaxIter = ' + str(self.maxiter) + ' exceeded. Task failed')
                task.status = 'F'
                continue

            if task.status in ('R', 'F', 'L'):
                # Skip running tasks
                if task.isrunning():
                    continue

                # Update and backup geometry and store cclib data if possible
                if not os.path.isfile(logfile):
                    print("No log file generated! Check " + task.full_gjf_path)
                    task.status = 'P'
                    continue

                last_ccdata, atoms = parse_gout(logfile)

                # Update energy and geometry info if at least one SCF is done
                if last_ccdata and hasattr(last_ccdata, 'scfenergies'):
                    if 'molstruct' in task.gjf:
                        task.old_coords = task.gjf['molstruct'].copy()
                    else:
                        task.old_coords = None
                    task.gjf['molstruct'] = atoms
                    task.ccdata = last_ccdata
                    last_ccdata.metadata['comments'] = []
                    for k in range(len(last_ccdata.atomcoords)):
                        last_ccdata.metadata['comments'].append('E_tot = {:6.5f}'.format(last_ccdata.scfenergies[k]))
                        last_ccdata.writexyz(add_index(task.folder + '/' + task.name + '.xyz', zerobased=True,
                                                       respect_file_extension=True), indices=k)
                    # Check if task is done
                    if checkflag(scenarios_dct['normal']['flags'], logfile, tail=scenarios_dct['normal']['tail']):
                        print(task.name + ': DONE')
                        task.status = 'D'
                        sh_execute('rm ' + task.folder + '/*.rwf')
                        continue
                elif hasattr(task, 'old_coords'):
                    task.gjf['molstruct'] = task.old_coords

                # Check if task is failed
                task_failed = False
                for fail_key, fail_val in scenarios_dct['fails'].items():
                    if checkflag(fail_val['flags'], logfile, tail=fail_val['tail']):
                        if dict2gauformat(task.gjf['route'], gau_route=True).lower() == "#p restart":
                            task.gjf = task.previousgjf.copy()
                        corrected_gjf = task.gjf.copy()
                        corrected_gjf.recurs_adjust(fail_val['corrector'])
                        if corrected_gjf == task.gjf and task.status in ['R', 'L']:
                            print(task.name + ': Task failed. Error: ' + fail_val['msg'])
                            task.status = 'F'
                        elif task.status == 'F':
                            print(task.name + ': Trying to revive failed task')
                            task.gjf = corrected_gjf
                            task.status = 'P'
                            task.k_iter += 1
                        continue

                # Prepare restart for erroneous tasks
                error_corrected = False
                for err_key, err_val in scenarios_dct['errors'].items():
                    if checkflag(err_val['flags'], logfile, tail=err_val['tail']):
                        if dict2gauformat(task.gjf['route'], gau_route=True).lower() == "#p restart":
                            task.gjf = task.previousgjf.copy()
                        print(task.name + ': Applying corrector for error:' + err_val['msg'])
                        task.gjf.recurs_adjust(err_val['corrector'])
                        task.status = 'P'
                        task.k_iter += 1
                        error_corrected = True
                if error_corrected:
                    task.previousgjf = task.gjf.copy()
                    continue

                # Two residual cases:
                if task.status in ['R', 'L']:
                    error_line = sh_execute("tail -30 " + logfile + " | grep 'Error termination via' ")
                    if error_line:  # 1. Unexpected error
                        print(task.name + ': Unexpected error: ' + error_line)
                        task.status = 'F'
                    else:  # 2. Simple restart
                        task.gjf.recurs_adjust(scenarios_dct['restart']['corrector'])
                        print(task.name + ': Restart')
                        task.status = 'P'

    def get_stats(self, stats_file='stats.txt', verb=False):

        stats = {'P': 0, 'R': 0, 'D': 0, 'F': 0, 'L':0}  # Pending, Running, Done, Failed, Loaded
        fullstory = []

        for task in self:
            stats[task.status] += 1
            task_story = '{:16s}'.format(task.name) + '\t' + '.' * task.k_iter + task.status
            fullstory.append(task_story)

        if verb:
            print('\n'.join(fullstory))

        if stats_file:
            with open(stats_file, 'a+') as stats_fid:
                for tsk in fullstory:
                    stats_fid.write(tsk + '\n')

        if stats['P'] == 0 and stats['R'] == 0 and stats['L'] == 0:
            print('No running or pending tasks!\n')
            if stats['F'] > 0:
                print('We\'ve done all we could with this input\n')
                print('To recover {} failed tasks, try to restart with altered scenarios'.format(stats['F']))
                sh_execute('touch STOP')
            else:
                print('Congratulations! All done!')
                sh_execute('touch DONE')

        print("Pending: {}, Running: {}, Done: {}, Failed: {}, Loaded: {}".format(
            stats['P'], stats['R'], stats['D'], stats['F'], stats['L']))
        return stats

    def dump(self, filename_pkl='database.pkl', filename_en='energies.txt'):
        with open(filename_en, 'w') as en_fid:
            en_fid.write('*' * 10 + '{:%Y-%m-%d %H:%M}'.format(datetime.now()) + '*' * 10 + '\n')
            en_fid.write('SCF energies in eV:\n')
            for task in self:
                if task.ccdata and hasattr(task.ccdata, 'scfenergies'):
                    en_fid.write('%-10s' % task.name + ': ' + '%-.10f' % (task.ccdata.scfenergies[-1]) + '\n')
                elif task.status == 'D':
                    last_ccdata, atoms = parse_gout(task.folder + '/' + task.out)
                    task.ccdata = last_ccdata
                    en_fid.write('%-10s' % task.name + ': ' + '%-.10f' % (task.ccdata.scfenergies[-1]) + '\n')
        with open(filename_pkl, 'wb') as pkl_fid:
            pickle.dump(self, pkl_fid)


# if __name__ == '__main__':
#     init_gjf = Gjf('../testfolder/complicated.gjf')
#
#     # jobtemplate = '../testfolder/job_arkuda.sh'
#     template = '../testfolder/job_oleg.sh'
#     dest_fold = '../DESTINATION'
#     poscars_file = '../testfolder/POSCARS'
#     database = GauCalcDB(scenarios='../scenarios.cjson', poscars_file=poscars_file, recalc_folder=dest_fold,
#                          machines_json='../machines.cjson', machine='rurik')
#     # database.submitjobs()
#     database.get_stats(verb=True)
#     database.dump()
#
#     del database
#
#     with open('database.pickle', 'rb') as db_fid:
#         db = pickle.load(db_fid)
#
#     # db
#
#     print(db)
#     # with open('pkl.pkl', 'rb') as f:
#     #     x = pickle.load(f)
#     # print(x[0].gjf)
#     # pass
