# import packages
from Visual import easy_time, visual_aa, visual_box, visual_mulbox
import itertools
import copy
from Extract import extract_row_plus, extract_col_plus, extract_trans


# AAC ##################################################################################################################
def feature_aac(pssm_aaid):
    aa = visual_aa()
    all_features = []
    for j in range(len(pssm_aaid)):
        line = pssm_aaid[j]
        aabox = visual_box(20)
        for i in line:
            if i in aa:
                aabox[aa.index(i)] += 1
        all_features.append(aabox)
    return all_features


# SAAC extract long sequence ###########################################################################################
def saac_el(eachfile):
    dn = 25
    dc = 10
    dn_1 = eachfile[0:4 * dn]
    dc_1 = eachfile[-dc:]
    dl_1 = eachfile[4 * dn:-dc]
    all_elfs = feature_aac([dn_1, dl_1, dc_1])
    out_box = []
    for i in all_elfs:
        out_box += i
    return out_box


# SAAC extract middle sequence
def saac_em(eachfile):
    dn = 25
    dc = 10
    dn_1 = eachfile[0:4 * dn]
    dc_1 = eachfile[-dc:]
    dl_1 = eachfile[-(dc + 20):-dc]
    all_elfs = feature_aac([dn_1, dl_1, dc_1])
    out_box = []
    for i in all_elfs:
        out_box += i
    return out_box


# SAAC extract short sequence
def saac_es(eachfile):
    dc = 10
    dn = int((len(eachfile) - dc) / 2)
    dn_1 = eachfile[0:dn]
    dc_1 = eachfile[-dc:]
    dl_1 = eachfile[dn:-dc]
    all_elfs = feature_aac([dn_1, dc_1, dl_1])
    out_box = []
    for i in all_elfs:
        out_box += i
    return out_box


# saac mian
def feature_saac(pssm_aaid):
    dn = 25
    dc = 10
    all_features = []
    for eachfile in pssm_aaid:
        each_fs = []
        if len(eachfile) >= (4 * dn + dc + 20):
            each_fs = saac_el(eachfile)
        if (4 * dn + dc) < len(eachfile) < (4 * dn + dc + 20):
            each_fs = saac_em(eachfile)
        if len(eachfile) <= (4 * dn + dc):
            each_fs = saac_es(eachfile)
        all_features.append(each_fs)
    return all_features


# Sliding window extract ###############################################################################################
def sw_extract(eachfile, eachaaid, lmda):
    eachmatrix = []
    sup_num = int((lmda - 1) / 2)
    sup_matrix = visual_mulbox(sup_num, len(eachfile[0]))
    sup_aaid = []
    for j in range(sup_num):
        sup_aaid.append('X')
    newfile = sup_matrix + eachfile + sup_matrix
    newaaid = sup_aaid + eachaaid + sup_aaid
    for j in range(sup_num, len(newfile) - sup_num):
        select_box = newfile[j - sup_num: j + sup_num + 1]
        eachmatrix.append([newaaid[j]] + select_box)
    return eachmatrix


# Sliding window plus
def sw_plus(eachmatrix):
    aa_index = visual_aa()
    matrix_400 = visual_mulbox(len(aa_index), len(eachmatrix[0][1]))
    for matrix in eachmatrix:
        if matrix[0] in aa_index:
            for line in matrix[1:]:
                for m in range(len(line)):
                    matrix_400[aa_index.index(matrix[0])][m] += line[m]
    for j in range(len(matrix_400)):
        for k in range(len(matrix_400[j])):
            matrix_400[j][k] = float('%.3f' % matrix_400[j][k])
    return matrix_400


# Sliding window reduce
def sw_reduce(raacode, eachplus):
    aa_index = visual_aa()
    reduce_fs = []
    for raa in raacode[1]:
        raa_box = raacode[0][raa]
        # 行合并
        score_box = visual_mulbox(len(raa_box), 20)
        for i in range(len(aa_index)):
            for j in range(len(raa_box)):
                if aa_index[i] in raa_box[j]:
                    score_box[j] = extract_row_plus(score_box[j], eachplus[i])
        # 列合并
        type_box = visual_mulbox(len(raa_box), len(raa_box))
        for i in range(len(aa_index)):
            for j in range(len(raa_box)):
                if aa_index[i] in raa_box[j]:
                    type_box[j] = extract_col_plus(type_box[j], score_box, i)
        type_box = extract_trans(type_box)
        out_box = []
        for i in type_box:
            out_box += i
        reduce_fs.append(out_box)
    return reduce_fs


# Sliding window main
def feature_sw(pssm_matrixes, pssm_aaid, raacode, lmda):
    sw_features = []
    start_e = 0
    for i in range(len(pssm_matrixes)):
        start_e += 1
        easy_time(start_e, len(pssm_matrixes))
        eachfile = pssm_matrixes[i]
        eachaaid = pssm_aaid[i]
        # 提取 L - lmda 个窗口
        eachmatrix = sw_extract(eachfile, eachaaid, int(lmda))
        # 合并20种aa为20x20维矩阵
        eachplus = sw_plus(eachmatrix)
        # 约化矩阵
        reduce_fs = sw_reduce(raacode, eachplus)
        sw_features.append(reduce_fs)
    return sw_features


# psekraac #############################################################################################################
def create_kmer_dict(simple_raa, kmer):
    kmer_dcit = {}
    for i in itertools.product("".join(simple_raa), repeat=kmer):
        j_str = ""
        for j in i:
            j_str = j_str + j
        kmer_dcit[j_str] = 0
    return kmer_dcit


# reduce
def kmer_reduce(raa_box, eachfile):
    out_box = []
    for i in eachfile:
        for j in raa_box:
            if i in j:
                out_box.append(j[0])
    simple_raa = []
    for i in raa_box:
        simple_raa.append(i[0])
    return [out_box, simple_raa]


def psekraac_reduce(raacode, raa, line):
    raa_box = raacode[0][raa]
    reduce_box = kmer_reduce(raa_box, line)
    reducefile = reduce_box[0]
    simple_raa = reduce_box[1]
    return reducefile, simple_raa


def kmer_count(simple_raa, line, kmer, gap, r):
    line = "".join(line)
    kmer_dcit = create_kmer_dict(simple_raa, kmer)
    kmer_dcit_c = kmer_dcit.copy()
    for i in range(0, len(line) - kmer*(r+1) + r + 1, gap + 1):
        kmer_dcit_c[line[i:i+kmer*(r+1):r+1]] = kmer_dcit_c[line[i:i+kmer*(r+1):r+1]] + 1
    return kmer_dcit_c


def kmer_change(data):
    out = []
    for key in data:
        out.append(data[key])
    return out


def kmer_part(raacode, pssm_aaid, k, g, lmda):
    sq_box = {}
    for i in range(len(pssm_aaid)):
        line = "".join(pssm_aaid[i])
        sq_box[i] = line
    out_box = []
    for line in pssm_aaid:
        mid_raa = []
        for raa in raacode[1]:
            reducefile, simple_raa = psekraac_reduce(raacode, raa, line)
            k_fs = kmer_count(simple_raa, reducefile, k, g, lmda)
            k_fs = kmer_change(k_fs)
            mid_raa.append(k_fs)
        out_box.append(mid_raa)
    return out_box


def psekraac_sum(k_fs, kg_fs, kgl_fs):
    out_box = []
    for i in range(len(k_fs)):
        raa_box_1 = k_fs[i]
        raa_box_2 = kg_fs[i]
        raa_box_3 = kgl_fs[i]
        new_box = []
        for j in range(len(raa_box_1)):
            dic_raa_1 = raa_box_1[j]
            dic_raa_2 = raa_box_2[j]
            dic_raa_3 = raa_box_3[j]
            new_dic = []
            for k in dic_raa_1:
                new_dic.append(k)
            for k in dic_raa_2:
                new_dic.append(k)
            for k in dic_raa_3:
                new_dic.append(k)
            new_box.append(new_dic)
        out_box.append(new_box)
    return out_box


# psekraac main
def feature_psekraac(raacode, pssm_aaid, k, g, lmda):
    # k
    k_fs = kmer_part(raacode, pssm_aaid, k, 0, 1)
    # k,g
    kg_fs = kmer_part(raacode, pssm_aaid, k, g, 1)
    # k,g,lmda
    kgl_fs = kmer_part(raacode, pssm_aaid, k, g, lmda)
    # sum
    psekraac_features = psekraac_sum(k_fs, kg_fs, kgl_fs)
    return psekraac_features


# kpssm DT #############################################################################################################
def kpssm_dt(eachfile, k):
    out_box = visual_box(len(eachfile[0]))
    for i in range(0, len(eachfile) - k):
        now_line = eachfile[i]
        next_line = eachfile[i + k]
        for j in range(len(now_line)):
            out_box[j] += now_line[j] * next_line[j]
    for i in range(len(out_box)):
        out_box[i] = out_box[i] / (len(eachfile) - k - 1)
    return out_box


# kpssm DDT
def kpssm_ddt(eachfile, k, aa_index):
    out_box = []
    for i in range((len(aa_index) - 1) * len(aa_index)):
        out_box.append(0)
    for i in range(0, len(eachfile) - k):
        now_line = eachfile[i]
        next_line = eachfile[i + k]
        n = -1
        for j in range(len(aa_index)):
            next_aa = copy.deepcopy(aa_index)
            now_aa = aa_index[j]
            next_aa.pop(next_aa.index(now_aa))
            for m in next_aa:
                n += 1
                out_box[n] += now_line[now_aa] * next_line[m]
    for i in range(len(out_box)):
        out_box[i] = out_box[i] / (len(eachfile) - k - 1)
    return out_box


# kpssm reduce
def kpssm_reduce(eachfile, raa_box):
    aa_index = visual_aa()
    # 列合并
    type_box = visual_mulbox(len(eachfile), len(raa_box))
    for li in range(len(eachfile)):
        for i in range(len(aa_index)):
            for j in range(len(raa_box)):
                if aa_index[i] in raa_box[j]:
                    type_box[li][j] += eachfile[li][i]
    # 列平均
    for i in range(len(type_box)):
        line = type_box[i]
        for j in range(len(line)):
            type_box[i][j] = type_box[i][j] / len(raa_box[j])
    return type_box


# k_gap_PSSM主程序
def feature_kpssm(pssm_matrixes, raacode):
    kpssm_features = []
    start_e = 0
    for eachfile in pssm_matrixes:
        start_e += 1
        easy_time(start_e, len(pssm_matrixes))
        mid_matrix = []
        for raa in raacode[1]:
            raa_box = raacode[0][raa]
            reducefile = kpssm_reduce(eachfile, raa_box)
            ddt_index = []
            for i in range(len(raa_box)):
                ddt_index.append(i)
            k = 3
            # DT
            dt_fs = kpssm_dt(reducefile, k)
            # DDT
            ddt_fs = kpssm_ddt(reducefile, k, ddt_index)
            mid_matrix.append(dt_fs + ddt_fs)
        kpssm_features.append(mid_matrix)
    return kpssm_features


# DT PSSM Reduce index #################################################################################################
def dtpssm_reduce(raa_box):
    out_box = []
    for i in raa_box:
        out_box.append(i[0])
    return out_box


# top 1 gram
def dtpssm_top_1(reducefile, reduce_index):
    out_box = []
    for line in reducefile:
        out_box.append(reduce_index[line.index(max(line))])
    return out_box


# d 0
def dtpssm_0(fs_top_1, reduce_index):
    out_box = visual_box(len(reduce_index))
    for i in fs_top_1:
        if i in reduce_index:
            out_box[reduce_index.index(i)] += 1
    return out_box


# d n
def dtpssm_n(fs_top_1, reduce_index, d):
    out_box = []
    out_index = []
    for i in reduce_index:
        for j in reduce_index:
            out_box.append(0)
            out_index.append(i + j)
    for i in range(len(fs_top_1) - d):
        if fs_top_1[i] + fs_top_1[i + d] in out_index:
            out_box[out_index.index(fs_top_1[i] + fs_top_1[i + d])] += 1
    return out_box


# PsePSSM主程序
def feature_dtpssm(pssm_matrixes, raacode, d):
    dt_features = []
    start_e = 0
    for eachfile in pssm_matrixes:
        start_e += 1
        easy_time(start_e, len(pssm_matrixes))
        mid_matrix = []
        for raa in raacode[1]:
            raa_box = raacode[0][raa]
            reducefile = kpssm_reduce(eachfile, raa_box)
            reduce_index = dtpssm_reduce(raa_box)
            # top_1_gram
            fs_top_1 = dtpssm_top_1(reducefile, reduce_index)
            # d_0
            fs_0 = dtpssm_0(fs_top_1, reduce_index)
            # d_n
            fs_n = []
            for m in range(1, d + 1):
                each_fs_n = dtpssm_n(fs_top_1, reduce_index, m + 1)
                fs_n += each_fs_n
            mid_matrix.append(fs_0 + fs_n)
        dt_features.append(mid_matrix)
    return dt_features
