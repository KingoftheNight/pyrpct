# import packages
import os
import math
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from Read import read_raac, read_feature, read_ssc, read_weblogo_main
from sklearn import svm
from sklearn import model_selection
from sklearn.metrics import roc_curve, auc
from Feature import create_kmer_dict
from pyecharts import options as opts
from pyecharts.charts import Bar, HeatMap, Pie, Sankey
from pyecharts.commons.utils import JsCode
from pyecharts.globals import ThemeType
from itertools import chain
from Visual import visual_box
plot_path = os.path.dirname(__file__)


# sort
def plot_eval_sort(box_v, box_i, type_b):
    out_v = []
    out_i = []
    for i in range(len(box_i)):
        out_v.append(box_v[i])
        out_i.append(int(box_i[i].strip(type_b)))
    n = len(out_i)
    for i in range(n):
        for j in range(0, n-i-1):
            if out_i[j] > out_i[j+1]:
                out_i[j], out_i[j+1] = out_i[j+1], out_i[j]
                out_v[j], out_v[j+1] = out_v[j+1], out_v[j]
    for i in range(len(out_i)):
        out_i[i] = type_b + str(out_i[i])
    return out_v, out_i


# 柱状图
def plot_eval_histogram(data, d_index, d_class, out, type_r):
    new_index = []
    for i in range(len(d_index)):
        new_index.append(d_index[i].strip(type_r))
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 汉显
    plt.rcParams['axes.unicode_minus'] = False  # 汉显
    plt.xlabel(d_class, fontsize=10)  # X轴标题
    plt.ylabel('ACC(%)', fontsize=10)  # Y轴标题
    plt.figure(figsize=(int(len(d_index)*0.2)+3, 5))
    plt.bar(new_index, data, width=0.8)  # 数据
    plt.title('ACC of each ' + d_class)  # 标题
    plt.grid(axis="y", c='g', linestyle='dotted')
    plt.savefig(os.path.join(out, 'ACC_' + d_class + '.png'), dpi=300)  # 保存
    plt.close()


# 密度图
def plot_eval_density(evaluate_score, out):
    values = []
    for key in evaluate_score:
        values.append(float(key[4]) * 100)
    s = pd.Series(values)
    sns.kdeplot(s, bins=10, hist=False, kde=True, axlabel='ACC')
    plt.savefig(os.path.join(out, 'ACC_Density.png'), dpi=300)
    plt.close()


# 热力图
def plot_eval_heatmap(evaluate_score, evaluate_key, out, t_s_index, s_s_index):
    data = {}
    for i in range(len(evaluate_score)):
        data[evaluate_key[i]] = float('%.4f' % (evaluate_score[i][4]))
    map_box = []
    min_num = float(data['t0s20']) * 100
    for key in data:
        if float(data[key]) * 100 > min_num:
            pass
        else:
            min_num = float(data[key]) * 100
    for s in s_s_index:
        mid_box = []
        for t in t_s_index:
            if t + s in data:
                mid_box.append(float('%.4f' % data[t + s]) * 100)
            else:
                mid_box.append(min_num - 10)
        map_box.append(mid_box)
    f, ax = plt.subplots(figsize=(int(len(t_s_index)*0.3)+4, 10))
    x = np.array(map_box)
    ax.set_title('ACC_Heatmap')
    ax.set_ylabel('Size')
    ax.set_xlabel('Type')
    sns.heatmap(x, cmap='YlGnBu', annot=True, mask=(x < min_num), vmax=100, linewidths=0.1, square=False,
                xticklabels=True, yticklabels=True)
    ax.set_xticklabels(t_s_index)
    ax.set_yticklabels(s_s_index)
    plt.savefig(os.path.join(out, 'ACC_Heatmap.png'), dpi=400)
    plt.close()


# evaluate plot
def plot_eval_main(cluster_t, t_index, out, cluster_s, s_index, evaluate_score, evaluate_key):
    # histogram
    cluster_s_t, t_s_index = plot_eval_sort(cluster_t, t_index, 't')
    plot_eval_histogram(cluster_s_t, t_s_index, 'Type', out, 't')
    cluster_s_s, s_s_index = plot_eval_sort(cluster_s, s_index, 's')
    plot_eval_histogram(cluster_s_s, s_s_index, 'Size', out, 's')
    # density
    plot_eval_density(evaluate_score, out)
    # heatmap
    plot_eval_heatmap(evaluate_score, evaluate_key, out, t_s_index, s_s_index)


# feature analize ######################################################################################################
def plot_fa_id(raacode, k):
    # raaPSSM
    raapssm = []
    for i in raacode:
        for j in raacode:
            raapssm.append(i[0] + j[0])
    # AAC
    aac = []
    aa = 'ARNDCQEGHILKMFPSTWYV'
    for i in aa:
        aac.append(i)
    # SAAC
    saac = []
    for i in range(3):
        for j in aa:
            saac.append(j)
    # raaKmer
    raa_kmer = []
    for i in raacode:
        raa_kmer.append(i[0])
    raakmer = list(create_kmer_dict(raa_kmer, k))
    # raaKPSSM
    raakpssm = []
    for i in raa_kmer:
        raakpssm.append(i)
    for i in raakmer:
        if i[0] != i[-1]:
            raakpssm.append(i)
    # raaSW
    raasw = []
    for i in raacode:
        for j in raacode:
            raasw.append(i[0] + j[0])
    # raaDTPSSM
    raadtpssm = []
    for i in raa_kmer:
        raadtpssm.append(i)
    for k in range(3):
        for i in raacode:
            for j in raacode:
                raadtpssm.append(i[0] + j[0])
    out_box = raapssm + aac + raakmer + raakpssm + saac + raasw + raadtpssm
    return out_box


def plot_fa_match(feature_id, fs_weight):
    out_dic = []
    n = 0
    for i in feature_id:
        i += i
        n += 1
        out_dic.append(fs_weight[str(n)])
    return out_dic


# 柱状图
def plot_fa_bar(data, all_data, index, out):
    c = (
        Bar()
        .add_xaxis(index)
        .add_yaxis("Amino Acids of AAC", data, itemstyle_opts=opts.ItemStyleOpts(color='gray'))
        .add_yaxis("All Amino Acids", all_data, itemstyle_opts=opts.ItemStyleOpts(color='orange'))
        .set_global_opts(
            title_opts=opts.TitleOpts(title="Amino Acids Features Weight"),
            yaxis_opts=opts.AxisOpts(name="weight"),
            xaxis_opts=opts.AxisOpts(name="type"),
            toolbox_opts=opts.ToolboxOpts(
                is_show=True, pos_top="top", pos_left="right", feature={"saveAsImage": {}})
        )
        .render(out)
    )
    return c


# 热图
def plot_fa_heatmap(x_index, y_index, value, out):
    c = (
        HeatMap(init_opts=opts.InitOpts(width="600px", height="600px"))
        .add_xaxis(x_index)
        .add_yaxis(
            "",
            y_index,
            value,
            label_opts=opts.LabelOpts(is_show=False, position="inside"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title="Features Weight of each pair of Amino Acids", pos_left='center'),
            visualmap_opts=opts.VisualMapOpts(is_show=False,
                                              max_=10, min_=0, range_color=['#F8F8FF', '#ADD8E6', '#000088']),
            toolbox_opts=opts.ToolboxOpts(
                is_show=True, pos_top="top", pos_left="right", feature={"saveAsImage": {}})
        )
        .render(out)
    )
    return c


# 饼图
def plot_fa_pie(data, out):
    c = (
        Pie()
        .add(
            "",
            data,
            radius=["30%", "75%"],
            center=["50%", "50%"],
            rosetype="radius",
            label_opts=opts.LabelOpts(is_show=True),
        )
        .set_global_opts(title_opts=opts.TitleOpts(title="The first 200-dimensional feature distribution",
                                                   pos_left='bottom'), toolbox_opts=opts.ToolboxOpts(
                is_show=True, pos_top="top", pos_left="right", feature={"saveAsImage": {}})
        )
        .render(out)
    )
    return c


def plot_fa_visual(feature_id, fw_index, fs_list, raacode, out):
    # 20 type
    aa = 'ARNDCQEGHILKMFPSTWYV'
    out_1 = {}
    for i in aa:
        out_1[i] = 0
    for i in feature_id:
        if len(i) == 1 and i in out_1:
            out_1[i] = out_1[i] + fw_index[feature_id.index(i)]
    aa20_data = []
    aa20_index = []
    for key in out_1:
        aa20_index.append(key)
        aa20_data.append(round(out_1[key], 0))
    # all type
    out_2 = {}
    for i in aa:
        out_2[i] = 0
    for i in feature_id:
        for j in i:
            if j in out_2:
                out_2[j] = out_2[j] + fw_index[feature_id.index(i)]
    allaa_data = []
    for key in out_2:
        allaa_data.append(round(out_2[key], 0))
    plot_fa_bar(aa20_data, allaa_data, aa20_index, os.path.join(out, 'Fw_bar.html'))
    # 20*20 type
    out_3 = {}
    for i in aa:
        for j in aa:
            out_3[i + j] = 0
    for i in feature_id:
        if i in out_3:
            out_3[i] = out_3[i] + fw_index[feature_id.index(i)]
    match_data = []
    match_index = []
    for i in aa:
        match_index.append(i)
        for j in aa:
            match_data.append([aa.index(i), aa.index(j), out_3[i+j]])
    plot_fa_heatmap(match_index, match_index, match_data, os.path.join(out, 'Fw_heatmap.html'))
    # 200 features
    le = len(raacode)
    data = [['raaPSSM', 0], ['raaKmer', 0], ['OAAC', 0], ['raaKPSSM', 0], ['SAAC', 0], ['raaSW', 0], ['raaDTPSSM', 0]]
    for i in range(200):
        if fs_list[i] < math.pow(le, 2):
            data[0][-1] = data[0][-1] + 1
        elif math.pow(le, 2) <= fs_list[i] < math.pow(le, 2) + 20:
            data[2][-1] = data[2][-1] + 1
        elif math.pow(le, 2) + 20 <= fs_list[i] < 2 * math.pow(le, 2) + 20:
            data[1][-1] = data[1][-1] + 1
        elif 2 * math.pow(le, 2) + 20 <= fs_list[i] < 3 * math.pow(le, 2) + 20:
            data[3][-1] = data[3][-1] + 1
        elif 3 * math.pow(le, 2) + 20 <= fs_list[i] < 3 * math.pow(le, 2) + 80:
            data[4][-1] = data[4][-1] + 1
        elif 3 * math.pow(le, 2) + 80 <= fs_list[i] < 4 * math.pow(le, 2) + 80:
            data[5][-1] = data[5][-1] + 1
        elif fs_list[i] >= 4 * math.pow(le, 2) + 80:
            data[-1][-1] = data[-1][-1] + 1
    plot_fa_pie(data, os.path.join(out, 'Fw_pie.html'))


def plot_fa_out(out, feature_id, fw_index):
    out_line = 'Feature ID\tFeature weight'
    for i in range(len(feature_id)):
        out_line += '\n' + feature_id[i] + '\t' + str(fw_index[i])
    with open(out, 'w') as f:
        f.write(out_line)


# fa main
def plot_fa_main(fs_weight, fs_list, raa, reduce, k, out):
    raa_file = os.path.join(os.path.join(plot_path, 'raacDB'), raa)
    raa_dict, raa_index = read_raac(raa_file)
    raacode = raa_dict[reduce]
    # feature id
    feature_id = plot_fa_id(raacode, k)
    # match
    fw_index = plot_fa_match(feature_id, fs_weight)
    # 特征可视化
    plot_fa_visual(feature_id, fw_index, fs_list, raacode, out)
    # 特征列表输出
    plot_fa_out(os.path.join(out, 'Feature_Weight.xls'), feature_id, fw_index)


# feature filter #######################################################################################################
# 特征分类
def filter_class(standard, size_fs, k):
    pssm_value = []
    aac_value = []
    kmer_value = []
    kpssm_value = []
    saac_value = []
    sw_value = []
    psepssm_value = []
    pmv = 0
    krv = 0
    acv = 0
    kmv = 0
    sav = 0
    swv = 0
    pev = 0
    for i in standard:
        if i < math.pow(size_fs, 2):
            pmv += 1
        if math.pow(size_fs, 2) <= i < (math.pow(size_fs, 2) + 20):
            acv += 1
        if (math.pow(size_fs, 2) + 20) <= i < (k * math.pow(size_fs, 2) + 20):
            krv += 1
        if (k * math.pow(size_fs, 2) + 20) <= i < ((k+1) * math.pow(size_fs, 2) + 20):
            kmv += 1
        if ((k+1) * math.pow(size_fs, 2) + 20) <= i < ((k+1) * math.pow(size_fs, 2) + 80):
            sav += 1
        if ((k+1) * math.pow(size_fs, 2) + 80) <= i < ((k+2) * math.pow(size_fs, 2) + 80):
            swv += 1
        if i >= ((k+2) * math.pow(size_fs, 2) + 80):
            pev += 1
        pssm_value.append(pmv / (1.5 * len(standard)))
        kmer_value.append(krv / (1.5 * len(standard)))
        aac_value.append(acv / (1.5 * len(standard)))
        kpssm_value.append(kmv / (1.5 * len(standard)))
        saac_value.append(sav / (1.5 * len(standard)))
        psepssm_value.append(pev / (1.5 * len(standard)))
        sw_value.append(swv / (1.5 * len(standard)))
    return pssm_value, aac_value, kmer_value, kpssm_value, saac_value, sw_value, psepssm_value


# 绘制折线图
def plot_filter_visual(data, pssm_value, aac_value, kmer_value, kpssm_value,
                       saac_value, sw_value, psepssm_value, type_p):
    x = []
    y = []
    for i in range(len(data)):
        x.append(i + 1)
    for j in data:
        y.append(float(j))
    plt.figure()
    plt.plot(x, y, label='ACC')
    plt.plot(x, pssm_value, color='blue', label='raaPSSM')
    plt.plot(x, kmer_value, color='green', label='raaKmer')
    plt.plot(x, aac_value, color='yellow', label='OAAC')
    plt.plot(x, kpssm_value, color='red', label='raaKPSSM')
    plt.plot(x, saac_value, color='pink', label='SAAC')
    plt.plot(x, sw_value, color='orange', label='raaSW')
    plt.plot(x, psepssm_value, color='gray', label='raaDTPSSM')
    plt.legend(bbox_to_anchor=(0., 1.09, 1., .102), loc=0, ncol=4, mode="expand", borderaxespad=0.)
    plt.xlabel("Feature Number")
    plt.ylabel("Acc")
    plt.title(type_p)
    max_x = y.index(max(y))
    max_y = max(y)
    max_pmv = pssm_value[max_x]
    max_krv = kmer_value[max_x]
    max_acv = aac_value[max_x]
    max_kmv = kpssm_value[max_x]
    max_sav = saac_value[max_x]
    max_swv = sw_value[max_x]
    max_pev = psepssm_value[max_x]
    plt.text(max_x, max_y, str(max_x + 1) + '(' + str(max_y * 100) + '%)', fontsize=10)
    plt.text(max_x, max_pmv, str(int(max_pmv * 1.5 * len(data))), fontsize=6)
    plt.text(max_x, max_krv, str(int(max_krv * 1.5 * len(data))), fontsize=6)
    plt.text(max_x, max_acv, str(int(max_acv * 1.5 * len(data))), fontsize=6)
    plt.text(max_x, max_kmv, str(int(max_kmv * 1.5 * len(data))), fontsize=6)
    plt.text(max_x, max_sav, str(int(max_sav * 1.5 * len(data))), fontsize=6)
    plt.text(max_x, max_swv, str(int(max_swv * 1.5 * len(data))), fontsize=6)
    plt.text(max_x, max_pev, str(int(max_pev * 1.5 * len(data))), fontsize=6)


# filter main
def plot_filter_main(in_path, relief_pool, fs_acc, k, out_path):
    size_fs = int(os.path.split(in_path)[-1].split('_')[0].split('s')[-1])
    # ifs
    pssm_value, aac_value, kmer_value, kpssm_value, saac_value, sw_value, psepssm_value = filter_class(relief_pool,
                                                                                                       size_fs, k)
    plot_filter_visual(fs_acc, pssm_value, aac_value, kmer_value, kpssm_value, saac_value,
                       sw_value, psepssm_value, 'IFS-Acc')
    plt.savefig(out_path, dpi=500, bbox_inches='tight')


# roc plot #############################################################################################################
def plot_roc_svm(af_data, af_label, c_number, ga):
    # 分割数据
    train_data, test_data, train_label, test_label = model_selection.train_test_split(af_data, af_label, test_size=.3,
                                                                                      random_state=0)
    # svm分类训练
    roc = svm.SVC(kernel='rbf', C=c_number, gamma=ga, probability=True)
    test_predict_label = roc.fit(train_data, train_label).decision_function(test_data)
    # roc坐标获取
    fpr, tpr, threshold = roc_curve(test_label, test_predict_label)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc


def plot_roc_line(fpr, tpr, roc_auc, out):
    plt.figure()
    plt.figure(figsize=(10, 10))
    lw = 2
    plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.savefig(out + '-roc.png', dpi=300)


# ROC main
def plot_roc_main(file, out, c_number, ga, now_path):
    np_data, np_label = read_feature(os.path.join(now_path, file))
    fpr, tpr, roc_auc = plot_roc_svm(np_data, np_label, float(c_number), float(ga))
    plot_roc_line(fpr, tpr, roc_auc, os.path.join(now_path, out))


# SSC ##################################################################################################################
# ssc cluster
def plot_ssc_cluster(source, target):
    sl, tl, vl = [], [], []
    for ti, taac in enumerate(target):
        taa_set = set(taac)
        aac_len = len(taac)
        for si, saac in enumerate(source):
            intersect = taa_set & set(saac)
            if intersect:
                sl.append(si)
                tl.append(ti)
                vl.append(len(intersect))
                aac_len -= len(intersect)
            if aac_len == 0:
                break
    return sl, tl, vl


# link
def plot_ssc_link(clusters):
    base_idx = 0
    source_idx, target_idx, values = [], [], []
    for i in range(len(clusters)-1):
        sl, tl, vl = plot_ssc_cluster(clusters[i], clusters[i+1])
        sidx = [i+base_idx for i in sl]
        base_idx += len(clusters[i])
        tidx = [i+base_idx for i in tl]
        source_idx.extend(sidx)
        target_idx.extend(tidx)
        values.extend(vl)
    return source_idx, target_idx, values


# sourse
def plot_ssc_sourse(labels, source_idx, target_idx, values):
    linkes = []
    for i in range(len(source_idx)):
        x_1 = source_idx[i]
        x_2 = target_idx[i]
        x_3 = values[i]
        mid_dic = {"source": labels[x_1], "target": labels[x_2], "value": x_3}
        if labels[x_1] != labels[x_2] and len(labels[x_1]) < len(labels[x_2]):
            linkes.append(mid_dic)
    return linkes


# nodes
def plot_ssc_nodes(linkes):
    name_box = []
    for dic in linkes:
        if dic["source"] not in name_box:
            name_box.append(dic["source"])
        if dic["target"] not in name_box:
            name_box.append(dic["target"])
    nodes = []
    for i in name_box:
        mid_dic = {"name": i}
        nodes.append(mid_dic)
    return nodes


# 绘制桑葚图
def plot_ssc_sankey(nodes, linkes, title_ssc, out):
    c = (
        Sankey()
        .add(
            title_ssc,
            nodes,
            linkes,
            linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source"),
            label_opts=opts.LabelOpts(position="right"),
        )
        .set_global_opts(title_opts=opts.TitleOpts(title="约化图谱"), toolbox_opts=opts.ToolboxOpts(
                is_show=True, pos_top="top", pos_left="right", feature={"saveAsImage": {}})
            )
        .render(out)
    )
    print("约化图谱保存于 " + c)


# ssc main
def plot_ssc_main(file, type_r, now_path):
    raa_file = os.path.join(os.path.join(plot_path, 'raacDB'), file)
    raac_list = read_ssc(raa_file, type_r)
    # get linkes
    source_idx, target_idx, values = plot_ssc_link(raac_list)
    labels = list(chain(*raac_list))
    linkes = plot_ssc_sourse(labels, source_idx, target_idx, values)
    # get nodes
    nodes = plot_ssc_nodes(linkes)
    # plot SSC
    title_ssc = "type" + type_r
    out = file + "_type" + type_r + "_SSC.html"
    out_path = os.path.join(now_path, out)
    plot_ssc_sankey(nodes, linkes, title_ssc, out_path)


# weblogo ##############################################################################################################
def plot_weblogo_check(line):
    add = 0
    for i in line:
        add += i
    if add > 100:
        mid = add - 100
        line[line.index(max(line))] = line[line.index(max(line))] - mid
    elif add < 100:
        mid = 100 - add
        line[line.index(max(line))] = line[line.index(max(line))] + mid
    return line


# weblogo change
def plot_weblogo_change(matrix, raacode):
    # 合并列
    new_matrix = []
    aa = 'ARNDCQEGHILKMFPSTWYV'
    for line in matrix:
        new_line = visual_box(len(raacode))
        for i in range(len(line)):
            for j in range(len(raacode)):
                if aa[i] in raacode[j]:
                    new_line[j] += line[i]
        new_line = plot_weblogo_check(new_line)
        new_matrix.append(new_line)
    # 提取列
    out_box = []
    for i in range(len(new_matrix[0])):
        mid_box = []
        for j in new_matrix:
            mid_box.append(j[i])
        out_box.append(mid_box)
    # 格式转换
    out_dic = []
    for each in out_box:
        mid_box = []
        for i in each:
            mid_dic = {"value": i, "percent": i / 100}
            mid_box.append(mid_dic)
        out_dic.append(mid_box)
    return out_dic


def plot_weblogo_draw(site_list, type_list, type_value, out):
    c = Bar(init_opts=opts.InitOpts(theme=ThemeType.LIGHT))
    c.add_xaxis(site_list)
    for i in range(len(type_list)):
        c.add_yaxis(type_list[i], type_value[i], stack="stack1", category_gap="50%")
    c.set_series_opts(
        label_opts=opts.LabelOpts(
            is_show=False,
            position="right",
            formatter=JsCode(
                "function(x){return Number(x.data.percent * 100).toFixed() + '%';}"
                ),
            )
        )
    if len(type_list) >= 8:
        c.set_global_opts(
            title_opts=opts.TitleOpts(title="Sequence Reduce Weblogo", pos_left="center"),
            datazoom_opts=[opts.DataZoomOpts(), opts.DataZoomOpts(type_="inside")],
            yaxis_opts=opts.AxisOpts(name="persernt(%)"),
            xaxis_opts=opts.AxisOpts(name="site"),
            toolbox_opts=opts.ToolboxOpts(is_show=True, pos_left="0px", pos_bottom="0px", feature={"saveAsImage": {}}),
            legend_opts=opts.LegendOpts(pos_right='0px', pos_top='10px', orient='vertical')
            )
    else:
        c.set_global_opts(
            title_opts=opts.TitleOpts(title="Sequence Reduce Weblogo", pos_left="center"),
            datazoom_opts=[opts.DataZoomOpts(), opts.DataZoomOpts(type_="inside")],
            yaxis_opts=opts.AxisOpts(name="persernt(%)"),
            xaxis_opts=opts.AxisOpts(name="site"),
            toolbox_opts=opts.ToolboxOpts(is_show=True, pos_left="0px", pos_bottom="0px", feature={"saveAsImage": {}}),
            legend_opts=opts.LegendOpts(pos_top='5%')
            )
    c.render(out)
    return c


# weblogo main
def plot_weblogo_main(file, raa, reduce, out, now_path):
    raa_file = os.path.join(os.path.join(plot_path, 'raacDB'), raa)
    raa_dict, raa_index = read_raac(raa_file)
    raacode = raa_dict[reduce]
    out_path = os.path.join(now_path, out + '.html')
    matrix = read_weblogo_main(file)
    type_value = plot_weblogo_change(matrix, raacode)
    site_list = []
    for i in range(len(type_value[0])):
        site_list.append(i+1)
    plot_weblogo_draw(site_list, raacode, type_value, out_path)
