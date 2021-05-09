# import packages
import os
import subprocess
import time
from Visual import visual_longcommand_linux, visual_longcommand_windows
blast_path = os.path.dirname(__file__)


# make database
def blast_makedb_linux(file, out):
    database_path = os.path.join(blast_path, 'blastDB')
    command = 'makeblastdb -in ' + file + ' -dbtype prot -parse_seqids -out ' + os.path.join(database_path, out)
    outcode = subprocess.Popen(command, shell=True)
    outcode.wait()


# 格式化数据库
def blast_makedb_windows(file, out):
    database_path = os.path.join(blast_path, 'blastDB')
    makedb_path = os.path.join(os.path.join(blast_path, 'bin'), 'makeblastdb.exe')
    command = makedb_path + ' -in ' + file + ' -dbtype prot -parse_seqids -out ' + os.path.join(database_path, out)
    outcode = subprocess.Popen(command, shell=True)
    outcode.wait()


# original psiblast
def blast_psiblast_linux(file, database, number, ev, out, now_path):
    file = os.path.join(os.path.join(now_path, 'Reads'), file)
    database_path = os.path.join(blast_path, 'blastDB')
    if 'PSSMs' not in os.listdir(now_path):
        root_pssm = os.path.join(now_path, 'PSSMs')
        os.makedirs(root_pssm)
    else:
        root_pssm = os.path.join(now_path, 'PSSMs')
    if out not in os.listdir(root_pssm):
        out_path = os.path.join(root_pssm, out)
        os.makedirs(out_path)
    else:
        out_path = os.path.join(root_pssm, out)
    order = 0
    for f in os.listdir(file):
        order += 1
        command = visual_longcommand_linux(file, f, database_path, database, number, ev, out_path)
        outcode = subprocess.Popen(command, shell=True)
        if outcode.wait() == 0:
            print('\r' + str(order) + '\tCompleted\t', end='', flush=True)
        else:
            print('\r' + str(order) + '\tProblems', end='', flush=True)
    if 'A' in os.listdir(now_path):
        os.remove('A')


# 调用psi-blast
def blast_psiblast_windows(file, database, number, ev, out, now_path):
    file = os.path.join(os.path.join(now_path, 'Reads'), file)
    database_path = os.path.join(blast_path, 'blastDB')
    psiblast_path = os.path.join(os.path.join(blast_path, 'bin'), 'psiblast.exe')
    if 'PSSMs' not in os.listdir(now_path):
        root_pssm = os.path.join(now_path, 'PSSMs')
        os.makedirs(root_pssm)
    else:
        root_pssm = os.path.join(now_path, 'PSSMs')
    if out not in os.listdir(root_pssm):
        path = os.path.join(root_pssm, out)
        os.makedirs(path)
    else:
        path = os.path.join(root_pssm, out)
    order = 0
    for f in os.listdir(file):
        order += 1
        command = visual_longcommand_windows(psiblast_path, file, f, database_path, database, number, ev, path)
        outcode = subprocess.Popen(command, shell=True)
        if outcode.wait() == 0:
            print('\r' + str(order) + '\tCompleted\t', end='', flush=True)
        else:
            print('\r' + str(order) + '\tProblems', end='', flush=True)
    if 'A' in os.listdir(now_path):
        os.remove('A')


# ray psiblast
def blast_rayblast_linux(folder, out, now_path):
    file_box = []
    out_box = []
    file_mid = []
    out_mid = []
    t = 0
    for f in os.listdir(os.path.join(os.path.join(now_path, 'Reads'), folder)):
        tp_path = os.path.join(os.path.join(now_path, 'Reads'), folder)
        pp_path = os.path.join(os.path.join(now_path, 'PSSMs'), out)
        if out not in os.listdir(os.path.join(now_path, 'PSSMs')):
            os.makedirs(pp_path)
        t += 1
        if t <= 20:
            file_mid.append(os.path.join(tp_path, f))
            out_mid.append(os.path.join(pp_path, f.split('.')[0]))
        else:
            file_box.append(file_mid)
            out_box.append(out_mid)
            file_mid = []
            out_mid = []
            t = 0
    start = time.time()
    for k in range(len(file_box)):
        out_file = ''
        for m in range(len(file_box[k])):
            out_file += file_box[k][m] + '@@' + out_box[k][m] + '\n'
        file_name = os.path.join(now_path, folder) + '_' + str(k)
        if folder + '_' + str(k) not in os.listdir(now_path):
            print(file_name)
            with open(file_name, 'w', encoding='UTF-8') as f:
                f.write(out_file)
            command = 'python Ray_blast.py ' + file_name
            subprocess.Popen(command, shell=True).communicate()
        else:
            pass
    if 'A' in os.listdir(now_path):
        os.remove('A')
    print("共计用时: {}s".format(time.time() - start))


def blast_rayblast_windows(folder, out, now_path):
    file_box = []
    out_box = []
    file_mid = []
    out_mid = []
    t = 0
    for f in os.listdir(os.path.join(os.path.join(now_path, 'Reads'), folder)):
        tp_path = os.path.join(os.path.join(now_path, 'Reads'), folder)
        pp_path = os.path.join(os.path.join(now_path, 'PSSMs'), out)
        if out not in os.listdir(os.path.join(now_path, 'PSSMs')):
            os.makedirs(pp_path)
        t += 1
        if t <= 20:
            file_mid.append(os.path.join(tp_path, f))
            out_mid.append(os.path.join(pp_path, f.split('.')[0]))
        else:
            file_box.append(file_mid)
            out_box.append(out_mid)
            file_mid = []
            out_mid = []
            t = 0
    file_box.append(file_mid)
    out_box.append(out_mid)
    start = time.time()
    for k in range(len(file_box)):
        out_file = ''
        for m in range(len(file_box[k])):
            out_file += file_box[k][m] + '@@' + out_box[k][m] + '\n'
        file_name = os.path.join(now_path, folder) + '_' + str(k)
        if folder + '_' + str(k) not in os.listdir(now_path):
            print(file_name)
            with open(file_name, 'w', encoding='UTF-8') as f:
                f.write(out_file)
            ray_command = os.path.join('pyrpct', 'Ray_blast_win.py')
            command = 'python ' + ray_command + ' ' + file_name
            subprocess.Popen(command, shell=True).communicate()
        else:
            pass
    if 'A' in os.listdir(now_path):
        os.remove('A')
    print("共计用时: {}s".format(time.time() - start))


# ray supplement
def blast_raysup_linux(folder, out, now_path):
    file_box = []
    out_box = []
    for f in os.listdir(os.path.join(os.path.join(now_path, 'Reads'), folder)):
        tp_path = os.path.join(os.path.join(now_path, 'Reads'), folder)
        pp_path = os.path.join(os.path.join(now_path, 'PSSMs'), out)
        if out not in os.listdir(os.path.join(now_path, 'PSSMs')):
            os.makedirs(pp_path)
        if f.split('.')[0] not in os.listdir(pp_path):
            file_box.append(os.path.join(tp_path, f))
            out_box.append(os.path.join(pp_path, f.split('.')[0]))
    start = time.time()
    out_file = ''
    for m in range(len(file_box)):
        out_file += file_box[m] + '@@' + out_box[m] + '\n'
    file_name = os.path.join(now_path, folder) + '_sup'
    if folder + '_sup' not in os.listdir(now_path):
        print(file_name)
        with open(file_name, 'w', encoding='UTF-8') as f:
            f.write(out_file)
        ray_command = os.path.join('pyrpct', 'Ray_blast.py')
        command = 'python ' + ray_command + ' ' + file_name
        subprocess.Popen(command, shell=True).communicate()
    else:
        pass
    if 'A' in os.listdir(now_path):
        os.remove('A')
    print("共计用时: {}s".format(time.time() - start))


def blast_raysup_windows(folder, out, now_path):
    file_box = []
    out_box = []
    for f in os.listdir(os.path.join(os.path.join(now_path, 'Reads'), folder)):
        tp_path = os.path.join(os.path.join(now_path, 'Reads'), folder)
        pp_path = os.path.join(os.path.join(now_path, 'PSSMs'), out)
        if out not in os.listdir(os.path.join(now_path, 'PSSMs')):
            os.makedirs(pp_path)
        if f.split('.')[0] not in os.listdir(pp_path):
            file_box.append(os.path.join(tp_path, f))
            out_box.append(os.path.join(pp_path, f.split('.')[0]))
    start = time.time()
    out_file = ''
    for m in range(len(file_box)):
        out_file += file_box[m] + '@@' + out_box[m] + '\n'
    file_name = os.path.join(now_path, folder) + '_sup'
    if folder + '_sup' not in os.listdir(now_path):
        print(file_name)
        with open(file_name, 'w', encoding='UTF-8') as f:
            f.write(out_file)
        ray_command = os.path.join('pyrpct', 'Ray_blast_win.py')
        command = 'python ' + ray_command + ' ' + file_name
        subprocess.Popen(command, shell=True).communicate()
    else:
        pass
    if 'A' in os.listdir(now_path):
        os.remove('A')
    print("共计用时: {}s".format(time.time() - start))
