# import packages
import os
import math
from Feature import *
from Visual import easy_time
from Read import read_raac, read_pssm
extract_path = os.path.dirname(__file__)


# 矩阵行相加
def extract_row_plus(data1, data2):
    new_data = []
    for i in range(len(data1)):
        new_data.append(data1[i] + data2[i])
    return new_data


# 矩阵列相加
def extract_col_plus(tp_box, score_box, i):
    for j in range(len(score_box)):
        line = score_box[j]
        tp_box[j] = tp_box[j] + line[i]
    return tp_box


# 转置
def extract_trans(type_box):
    new_box = []
    for i in range(len(type_box)):
        col = []
        for j in type_box:
            col.append(j[i])
        new_box.append(col)
    return new_box


# 获取raac
def extract_raa(raa, outfolder, now_path):
    # 获取氨基酸约化密码表
    raa_path = os.path.join(extract_path, 'raacDB')
    raa_file = os.path.join(raa_path, raa)
    if raa in os.listdir(raa_path):
        raacode = read_raac(raa_file)
        if outfolder not in os.listdir(now_path):
            outfolder = os.path.join(now_path, outfolder)
            os.makedirs(outfolder)
        else:
            outfolder = os.path.join(now_path, outfolder)
    else:
        with open(raa_file, 'w') as f:
            f.write('type 1 size ' + str(len(raa.split('-'))) + ' ' + raa)
        raacode = read_raac(raa_file)
        if outfolder not in os.listdir(now_path):
            outfolder = now_path
        else:
            outfolder = now_path
        os.remove(raa_file)
    return raacode, outfolder


# str to float
def extract_change(pssm_matrixes):
    out_box = []
    for i in range(len(pssm_matrixes)):
        mid_box = []
        for j in range(len(pssm_matrixes[i])):
            next_box = []
            for k in range(len(pssm_matrixes[i][j])):
                next_box.append(float(pssm_matrixes[i][j][k]))
            mid_box.append(next_box)
        out_box.append(mid_box)
    return out_box


# 提取矩阵特征
def extract_features(pssm_matrixes, pssm_aaid):
    all_features = []
    start_e = 0
    for i in range(len(pssm_matrixes)):
        start_e += 1
        easy_time(start_e, len(pssm_matrixes))
        each_matrix = pssm_matrixes[i]
        matrix_400 = []
        aa_index = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        for aa in aa_index:
            aa_score = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0]
            for j in range(len(each_matrix)):
                line = each_matrix[j]
                if pssm_aaid[i][j] == aa:
                    for k in range(len(line)):
                        aa_score[k] = aa_score[k] + line[k]
            matrix_400.append(aa_score)
        all_features.append(matrix_400)
    return all_features


# long to short
def extract_short(pssm_features):
    out_box = []
    for i in range(len(pssm_features)):
        mid_box = []
        for j in range(len(pssm_features[i])):
            next_box = []
            for k in range(len(pssm_features[i][j])):
                next_box.append(float('%.3f' % pssm_features[i][j][k]))
            mid_box.append(next_box)
        out_box.append(mid_box)
    return out_box


# 约化矩阵
def extract_reduce(pssm_features, raacode, pssm_type):
    all_features = []
    aa_index = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    start_e = 0
    for raa in raacode[1]:
        start_e += 1
        easy_time(start_e, len(raacode[1]))
        raa_box = raacode[0][raa]
        mid_box = []
        for k in range(len(pssm_features)):
            eachfile = pssm_features[k]
            eachtype = pssm_type[k]
            # 行合并
            score_box = visual_mulbox(len(raa_box), 20)
            for i in range(len(aa_index)):
                for j in range(len(raa_box)):
                    if aa_index[i] in raa_box[j]:
                        score_box[j] = extract_row_plus(score_box[j], eachfile[i])
            # 列合并
            type_box = visual_mulbox(len(raa_box), len(raa_box))
            for i in range(len(aa_index)):
                for j in range(len(raa_box)):
                    if aa_index[i] in raa_box[j]:
                        type_box[j] = extract_col_plus(type_box[j], score_box, i)
            # 转置
            type_box = extract_trans(type_box)
            mid_box.append([eachtype, type_box])
        all_features.append(mid_box)
    return all_features


# 矩阵转置
def extract_transform(data):
    new_data = []
    for i in data:
        new_data += i
    return new_data


# 归一化
def extract_scale(data):
    new_data = []
    for i in data:
        if i >= - 709:
            f = 1 / (1 + math.exp(-i))
            new_data.append(f)
        else:
            f = 1 / (1 + math.exp(-709))
            new_data.append(f)
        # new_data.append('%.6f' % f)
    return new_data


# 保存特征
def extract_save(raa_features, aac_features, psekraac_features, kpssm_features, saac_features, sw_features,
                 dtpssm_features, outfolder, raa_list):
    start_e = 0
    for k in range(len(raa_features)):
        start_e += 1
        easy_time(start_e, len(raa_features))
        eachraa = raa_features[k]
        out_file = ''
        for i in range(len(eachraa)):
            eachfile = eachraa[i]
            eachaac = aac_features[i]
            eachkpssm = kpssm_features[i][k]
            eachkmer = psekraac_features[i][k]
            eachsaac = saac_features[i]
            eachsw = sw_features[i][k]
            eachdtpssm = dtpssm_features[i][k]
            type_m = eachfile[0]
            data_m = eachfile[1]
            data_m = extract_transform(data_m) + eachaac + eachkmer + eachkpssm + eachsaac + eachsw + eachdtpssm
            # 归一化
            data_m = extract_scale(data_m)
            mid_file = type_m
            for j in range(len(data_m)):
                mid_file += ' ' + str(j + 1) + ':' + str(data_m[j])
            out_file += mid_file + '\n'
        path = os.path.join(outfolder, raa_list[k] + '_rpct.fs')
        with open(path, 'w') as f2:
            f2.write(out_file)
            f2.close()


# extract main
def extract_main(positive, negative, outfolder, raa, lmda, now_path):
    # 获取raac
    raacode, outfolder = extract_raa(raa, outfolder, now_path)
    # 处理地址
    positive = os.path.join(os.path.join(now_path, 'PSSMs'), positive)
    negative = os.path.join(os.path.join(now_path, 'PSSMs'), negative)
    # positive
    pssm_matrixes, pssm_aaid, pssm_type = read_pssm(positive, [], [], [], '0')
    # negative
    pssm_matrixes, pssm_aaid, pssm_type = read_pssm(negative, pssm_matrixes, pssm_aaid, pssm_type, '1')
    # PSSM特征提取
    pssm_matrixes = extract_change(pssm_matrixes)
    pssm_features = extract_features(pssm_matrixes, pssm_aaid)
    pssm_features = extract_short(pssm_features)
    # Sequence PseKRAAC
    k, g, m = 2, 1, 1
    psekraac_features = feature_psekraac(raacode, pssm_aaid, k, g, m)
    # AAC特征提取
    aac_features = feature_aac(pssm_aaid)
    # SAAC特征提取
    saac_features = feature_saac(pssm_aaid)
    # kpssm特征提取
    kpssm_features = feature_kpssm(pssm_matrixes, raacode)
    # psepssm特征提取
    dtpssm_features = feature_dtpssm(pssm_matrixes, raacode, 3)
    # PSSM滑窗
    sw_features = feature_sw(pssm_matrixes, pssm_aaid, raacode, lmda)
    # 矩阵约化
    raa_features = extract_reduce(pssm_features, raacode, pssm_type)
    # 生成特征文件
    extract_save(raa_features, aac_features, psekraac_features, kpssm_features, saac_features, sw_features,
                 dtpssm_features,
                 outfolder, raacode[1])
