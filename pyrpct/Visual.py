# import packages
import os


# time
def easy_time(start_e, end_e):
    print('\r>>>' + str(start_e) + "~" + str(end_e), end='', flush=True)


# 20 aa
def visual_aa():
    return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


# single level box
def visual_box(i):
    out = []
    for k in range(i):
        out.append(0)
    return out


# multy level box
def visual_mulbox(i, j):
    out = []
    for k in range(i):
        mid = []
        for m in range(j):
            mid.append(0)
        out.append(mid)
    return out


# block
def visual_block(any_thing):
    any_thing += any_thing


# long blast command
def visual_longcommand_linux(file, f, database_path, database, number, ev, out_path):
    command = 'psiblast -query ' + os.path.join(
        file, f) + ' -db ' + os.path.join(
        database_path, database
    ) + ' -num_iterations ' + number + ' -evalue ' + ev + ' -out A' + ' -out_ascii_pssm ' + os.path.join(
        out_path, f.split('.')[0])
    return command


def visual_longcommand_windows(psiblast_path, file, f, database_path, database, number, ev, path):
    command = psiblast_path + ' -query ' + os.path.join(
        file, f) + ' -db ' + os.path.join(
        database_path, database
    ) + ' -num_iterations ' + number + ' -evalue ' + ev + ' -out A' + ' -out_ascii_pssm ' + os.path.join(
        path, f.split('.')[0])
    return command


def visual_longcommand_raylinux(file, database_path, out):
    command = 'psiblast -query ' + file + ' -db ' + database_path + ' -num_iterations 3 -evalue 0.001 -out A -out_ascii_pssm ' + out
    return command


def visual_longcommand_raywindows(psiblast_path, file, database_path, out):
    command = psiblast_path + ' -query ' + file + ' -db ' + database_path + ' -num_iterations 3 -evalue 0.001 -out A -out_ascii_pssm ' + out
    return command
