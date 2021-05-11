# import packages
import os
import math
import random
from Plot import plot_fa_main, plot_filter_main
from SVM import svm_evaluate
from TES import tes_score
from Read import read_pca
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from Visual import easy_time
import numpy as np


# ojld distance
def ojld_distance(p_data, n_data, test_data, number):
    p_distance = []
    n_distance = []
    for line in p_data:
        p_distance.append(math.pow((test_data[number] - line[number]), 2))
    for line in n_data:
        n_distance.append(math.pow((test_data[number] - line[number]), 2))
    return [min(p_distance), min(n_distance)]


# relief method
def filter_relief(number, feature_class, feature_line, cycle):
    feature_standard = []
    for line in feature_line:
        mid_box = []
        for i in line:
            mid_box.append(float(i.split(':')[-1]))
        feature_standard.append(mid_box)
    p_data = []
    n_data = []
    for i in range(len(feature_standard)):
        if feature_class[i] == '0':
            p_data.append(feature_standard[i])
        if feature_class[i] == '1':
            n_data.append(feature_standard[i])
    weight = 0
    m = 0
    for m in range(cycle):
        rand_num = random.randint(0, len(feature_standard) - 1)
        if feature_class[rand_num] == '0':
            distance_box = ojld_distance(p_data, n_data, feature_standard[rand_num], number)
            weight += -distance_box[0] + distance_box[1]
        if feature_class[rand_num] == '1':
            distance_box = ojld_distance(p_data, n_data, feature_standard[rand_num], number)
            weight += -distance_box[1] + distance_box[0]
    aver_weight = weight / (m + 1)
    aver_weight = 1 / (1 + math.exp(-aver_weight))
    return aver_weight


# 概率计算
def filter_fscore(number, feature_class, feature_line):
    type_both = 0
    type_a = 0
    type_b = 0
    t0 = 0
    t1 = 0
    for i in range(len(feature_class)):
        if feature_class[i] == '0':
            type_a += float(feature_line[i][number].split(':')[-1])
            t0 += 1
        else:
            type_b += float(feature_line[i][number].split(':')[-1])
            t1 += 1
        type_both += float(feature_line[i][number].split(':')[-1])
    avg_0 = type_a / t0
    avg_1 = type_b / t1
    avg_both = type_both / len(feature_class)
    f_son = math.pow(avg_0 - avg_both, 2) + math.pow(avg_1 - avg_both, 2)
    avg_m_0 = 0
    avg_m_1 = 0
    for i in range(len(feature_class)):
        if feature_class[i] == '0':
            avg_m_0 += (math.pow(float(feature_line[i][number].split(':')[-1]) - avg_0, 2))
        else:
            avg_m_1 += (math.pow(float(feature_line[i][number].split(':')[-1]) - avg_1, 2))
    f_mother = avg_m_0 / (t0 - 1) + avg_m_1 / (t1 - 1)
    if f_mother != 0:
        f_score = f_son / f_mother
    else:
        f_score = -0.1
    return f_score


# 排序
def filter_sort(data):
    arr = []
    for i in data:
        arr.append(i)
    index = []
    for i in range(len(arr)):
        index.append(i)
    for i in range(len(arr) - 1):
        min_index = i
        for j in range(i + 1, len(arr)):
            if arr[j] < arr[min_index]:
                min_index = j
        index[min_index], index[i] = index[i], index[min_index]
        arr[min_index], arr[i] = arr[i], arr[min_index]
    # 倒序输出
    re_index = []
    for i in range(len(index) - 1, -1, -1):
        re_index.append(index[i])
    return re_index


# 特征组合测试
def filter_test(d, feature, c_number, gamma, crossv, now_path):
    fs_acc = []
    filter_data = []
    for k in d:
        filter_data.append(k[0])
    start_e = 0
    for i in range(len(feature)):
        start_e += 1
        key = feature[i]
        for j in range(len(d)):
            filter_data[j] += ' ' + str(i + 1) + ':' + d[j].split(' ')[key + 1].strip('\n').split(':')[-1]
        out_content = ''
        for n in filter_data:
            out_content += n + '\n'
        with open('mid-ifs', 'w') as ot:
            ot.write(out_content)
        test_label, predict_label = svm_evaluate(os.path.join(now_path, 'mid-ifs'),
                                                 float(c_number), float(gamma), int(crossv))
        standard_num = tes_score(test_label, predict_label)
        single_acc = str('%.4f' % (standard_num[4]))
        fs_acc.append(single_acc)
        os.remove('./mid-ifs')
        easy_time(start_e, len(feature))
    return fs_acc


# filter main
def filter_main(in_path, out_path, c_number, gamma, crossv, cycle, raa, reduce, now_path):
    in_path = os.path.join(now_path, in_path)
    if out_path not in os.listdir(now_path):
        out_path = os.path.join(now_path, out_path)
        os.makedirs(out_path)
    else:
        out_path = os.path.join(now_path, out_path)
    with open(in_path, 'r') as in_file:  # 读取文件
        d = in_file.readlines()
        in_file.close()
    feature_class = []
    feature_line = []
    for i in range(len(d)):
        lines = d[i].strip('\n')
        feature_class.append(lines.split(' ')[0])
        mid_box = lines.split(' ')[1:]
        if len(mid_box[-1]) == 0:
            mid_box = mid_box[:-1]
        # mid_box = filter_check(mid_box, num)
        feature_line.append(mid_box)
        out = lines.split(' ')[0] + ' '
        for j in mid_box:
            out += j + ' '
        d[i] = out + '\n'
    relief_list = []
    start_num = 0
    for each_number in range(len(feature_line[0])):
        start_num += 1
        easy_time(start_num, len(feature_line[0]))
        type_relief = filter_relief(each_number, feature_class, feature_line, int(cycle))  # 求得每个特征relief
        type_fscore = filter_fscore(each_number, feature_class, feature_line)  # 求得每个特征f-score
        complex_num = 1 / (math.exp(-type_relief) + 1) + 1 / (math.exp(-type_fscore) + 1)
        relief_list.append(complex_num)
    relief_pool = filter_sort(relief_list)  # 排序
    # fsw画图
    k = 2  # 测试
    out_dic = {}
    for i in relief_pool:
        out_dic[str(i + 1)] = relief_list[i]
    fs_weight, fs_list = out_dic, relief_pool
    plot_fa_main(fs_weight, fs_list, raa, reduce, k, out_path)
    # 折线图
    fs_acc = filter_test(d, relief_pool, c_number, gamma, crossv, now_path)
    print('\n特征筛选完成，导出结果中...')
    plot_filter_main(in_path, relief_pool, fs_acc, k, os.path.join(out_path, 'ACC-ifs.png'))
    ifs_content = 'IFS-feature-sort: '
    for i in relief_pool:
        ifs_content += str(i + 1) + ' '
    with open(os.path.join(out_path, 'Fsort-ifs.txt'), 'w') as f4:
        f4.write(ifs_content)
        f4.close()


# PCA method ###########################################################################################################
def pca_scale(x_train):
    min_max_scaler = MinMaxScaler(feature_range=(-1, 1))
    scaler = min_max_scaler.fit(x_train)
    x_train_ = scaler.transform(x_train)
    return x_train_


# feature selection
def pca_selection(data):
    x = data.iloc[:, 1:]
    x_s = pca_scale(x)
    pca = PCA(n_components=1)
    pca.fit(x_s)
    pc1_loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    pc1_featurescore = pd.DataFrame({'Feature': x.columns,
                                     'PC1_loading': pc1_loadings.T[0],
                                     'PC1_loading_abs': abs(pc1_loadings.T[0])})
    pc1_featurescore = pc1_featurescore.sort_values('PC1_loading_abs', ascending=False)
    feature_selection = []
    for i in pc1_featurescore['Feature']:
        feature_selection.append(i - 1)
    return feature_selection


# PCA main
def pca_main(in_path, out_path, c_number, gamma, crossv, now_path):
    # 读取文件
    in_path = os.path.join(now_path, in_path)
    if out_path not in os.listdir(now_path):
        out_path = os.path.join(now_path, out_path)
        os.makedirs(out_path)
    else:
        out_path = os.path.join(now_path, out_path)
    features_data, data = read_pca(in_path)
    # PCA
    features_array = pd.DataFrame(features_data)
    fs_matrix = pca_selection(features_array)
    # 调用filter
    fs_acc = filter_test(data, fs_matrix, c_number, gamma, crossv, now_path)
    print('\n特征筛选完成，导出结果中...')
    # 绘图
    k = 2  # 测试
    plot_filter_main(in_path, fs_matrix, fs_acc, k, os.path.join(out_path, 'Fsort-pca'))
    out_file = 'IFS-feature-sort: '
    for j in fs_matrix:
        out_file += str(j + 1) + ' '
    with open(os.path.join(out_path, 'Fsort-pca.txt'), 'w', encoding='UTF-8') as f:
        f.write(out_file)
        f.close()
