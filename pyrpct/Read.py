# import packages
import os
import copy
import numpy as np
import sys
read_path = os.path.dirname(__file__)
sys.path.append(read_path)
from Visual import easy_time, visual_block, visual_box


# functions
def read_fasta(file, out, now_path):
    if 'Reads' not in os.listdir(now_path):
        root_read = os.path.join(now_path, 'Reads')
        os.makedirs(root_read)
    else:
        root_read = os.path.join(now_path, 'Reads')
    with open(file, 'r') as u:
        lines = u.readlines()
    result = ''
    for i in lines:
        i = i.strip()
        if i:
            if i[0] == '>':
                result = result + '\n' + i + '\n'
            else:
                result = result + i
    result = result[1:]
    order = 0
    filelist = result.split('\n')
    if out not in os.listdir(root_read):
        out_path = os.path.join(root_read, out)
        os.makedirs(out_path)
    else:
        out_path = os.path.join(root_read, out)
    for line in range(len(filelist)):
        if '>' in filelist[line]:
            order += 1
            if line < len(filelist):
                with open(os.path.join(out_path, str(order) + '.fasta'), 'w') as f:
                    f.write(filelist[line] + '\n' + filelist[line + 1])


def read_raac(file):
    with open(file, 'r') as code:
        raacode = code.readlines()
    raa_dict = {}
    raa_index = []
    for eachline in raacode:
        each_com = eachline.strip('\n').split()
        raa_com = each_com[-1].split('-')
        raa_ts = 't' + each_com[1] + 's' + each_com[3]
        raa_dict[raa_ts] = raa_com
        raa_index.append(raa_ts)
    return raa_dict, raa_index


def read_database(file_name, now_path):
    f = open(os.path.join(now_path, file_name), 'r', encoding='UTF-8')
    namebox = []
    eachsequence = ''
    line = 1
    t = 0
    while line:
        t += 1
        line = f.readline()
        if '>' in line:
            if line.split(' ')[0] not in namebox and line.split(' ')[0].upper() not in namebox:
                namebox.append(line.split(' ')[0].upper())
                eachsequence += line.split(' ')[0].upper() + '\n'
                writeable = 'Ture'
            else:
                writeable = 'False'
        else:
            if writeable == 'Ture':
                eachsequence += line
            else:
                pass
        if t == 20000:
            with open(os.path.join(now_path, 'ND_' + file_name), 'a', encoding='UTF-8') as o:
                o.write(eachsequence)
                o.close()
            eachsequence = ''
            t = 0
    with open(os.path.join(now_path, 'ND_' + file_name), 'a', encoding='UTF-8') as o:
        o.write(eachsequence)
        o.close()
    f.close()


def read_pssm(path, pssm_matrixes, pssm_aaid, pssm_type, type_p):
    start_e = 0
    for i in range(len(os.listdir(path))):
        start_e += 1
        easy_time(start_e, len(os.listdir(path)))
        eachfile = os.listdir(path)[i]
        with open(os.path.join(path, eachfile), 'r') as f1:
            data = f1.readlines()
        matrix = []
        aa_id = []
        end_matrix = 0
        for j in data:
            if 'Lambda' in j and 'K' in j:
                end_matrix = data.index(j)
                break
        for eachline in data[3:end_matrix - 1]:
            row = eachline.split()
            newrow = row[0:22]
            for k in range(2, len(newrow)):
                newrow[k] = int(newrow[k])
            nextrow = newrow[2:]
            matrix.append(nextrow)
            aa_id.append(newrow[1])
        pssm_matrixes.append(matrix)
        pssm_aaid.append(aa_id)
        pssm_type.append(type_p)
    return copy.deepcopy(pssm_matrixes), copy.deepcopy(pssm_aaid), copy.deepcopy(pssm_type)


def read_filter(file, filter_index, out, number, now_path):
    # 读取特征排序以及特征文件
    with open(filter_index, 'r', encoding='UTF-8') as f1:
        data = f1.readlines()
    index = data[0].split(' ')[1:-1]
    with open(file, 'r', encoding='UTF-8') as f2:
        data = f2.readlines()
    # 提取矩阵特征
    type_f = []
    matrix = []
    for line in data:
        line = line.split(' ')
        type_f.append(line[0])
        mid_box = line[1:]
        for i in range(len(mid_box)):
            mid_box[i] = mid_box[i].split(':')[-1]
        matrix.append(mid_box)
    # 生成特征筛选文件
    new_matrix = ''
    for i in range(len(matrix)):
        line = matrix[i]
        mid_file = type_f[i]
        order = 0
        for j in range(0, int(number)):
            key = index[j]
            order += 1
            mid_file += ' ' + str(order) + ':' + line[int(key) - 1].strip('\n')
        new_matrix += mid_file + '\n'
    with open(os.path.join(now_path, out + '-fffs_rpct.fs'), 'w') as f3:
        f3.write(new_matrix[:-1])


def read_feature(file_name):
    with open(file_name, 'r') as f1:
        file = f1.readlines()
    # 提取特征list
    features = []
    features_label = []
    for i in file:
        line = i.strip('\n').split(' ')
        fs_box = line[1:]
        mid_box = []
        for j in fs_box:
            mid_box.append(float(j.split(':')[-1]))
        features.append(mid_box)
        features_label.append(int(line[0]))
    # 转换为数组
    np_data = np.array(features)
    np_label = np.array(features_label)
    return np_data, np_label


def read_hys(file):
    with open(file, 'r') as f:
        data = f.readlines()
    cg_box = {}
    for i in data:
        cg_box[i.split('\t')[0]] = [i.split('\t')[1].split(': ')[-1], i.split('\t')[2].split(': ')[-1]]
    return cg_box


def read_intlen_hys(fs_key, cg_file):
    with open(cg_file, 'r') as hf:
        data = hf.readlines()
        hf.close()
    cg_box = visual_box(len(fs_key))
    for each in data:
        eachline = each.split('	')
        c_num = float(eachline[1].split(': ')[-1])
        gamma = float(eachline[2].strip('\n').split(': ')[-1])
        if eachline[0] in fs_key:
            cg_box[fs_key.index(eachline[0])] = [c_num, gamma]
    return cg_box


def read_intro(name):
    with open(os.path.join(os.path.join(read_path, 'bin'), 'intro'), 'r', encoding='UTF-8') as f:
        data = f.readlines()
    out_line = ''
    writable = 'no'
    for i in range(len(data)):
        line = data[i].strip('\n')
        if '>' in line:
            if name in line:
                writable = 'yes'
            else:
                writable = 'no'
        if writable == 'yes':
            out_line += data[i].strip('\n') + '\n'
    return out_line[len(name)+2:-1]


def read_pca(in_path):
    with open(in_path, 'r') as f1:
        data = f1.readlines()
        f1.close()
    out_box = []
    for line in data:
        each_box = line.strip('\n').split(' ')
        mid_box = [float(each_box[0])]
        for i in each_box[1:]:
            mid_box.append(float(i.split(':')[-1]))
        out_box.append(mid_box)
    out_box = np.array(out_box)
    return out_box, data


def read_ssc(raa_file, type_r):
    with open(raa_file, "r") as f:
        data = f.readlines()
    out_box = []
    for line in data:
        line = line.strip("\n").split(" ")
        if line[1] == type_r:
            out_box.append(line[4])
    all_sq = ""
    for i in out_box[0]:
        if i != "-":
            all_sq += i + "-"
    out_box.append(all_sq[:-1])
    for i in range(len(out_box)):
        out_box[i] = out_box[i].split("-")
    return out_box[::-1]


def read_weblogo_score(newrow):
    out_box = []
    a = 0
    for i in newrow:
        i = int(i)
        if i > 0:
            a += i
    for i in newrow:
        i = int(i)
        if i > 0:
            out_box.append(i * 100 / a)
        else:
            out_box.append(0)
    return out_box


def read_weblogo_main(file):
    with open(file, 'r') as f:
        data = f.readlines()
    end_matrix = 0
    for j in data:
        if 'Lambda' in j and 'K' in j:
            end_matrix = data.index(j)
            break
    matrix = []
    for eachline in data[3:end_matrix - 1]:
        row = eachline.split()
        newrow = row[2:22]
        newrow = read_weblogo_score(newrow)
        for i in range(len(newrow)):
            newrow[i] = int(newrow[i])
        matrix.append(newrow)
    return matrix


def read_precaution():
    with open(os.path.join(read_path, 'README'), 'r', encoding='UTF-8') as f:
        prec_data = f.read()
    return prec_data
