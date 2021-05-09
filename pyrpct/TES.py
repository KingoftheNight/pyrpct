# import packages
import os
from SVM import svm_train, svm_evaluate, svm_grid, svm_predict, svm_acc
from Read import read_hys, read_intro
from Visual import easy_time, visual_box, visual_block
from Plot import plot_eval_main
from math import sqrt
tes_path = os.path.dirname(__file__)


# make Hyperparameters
def tes_makehys(folder, c, g, out, now_path):
    out_file = ''
    folder = os.path.join(now_path, folder)
    out = os.path.join(now_path, out)
    for i in os.listdir(folder):
        out_file += i + '	C_numbr: ' + c + '	Gamma: ' + g + '\n'
    with open(out, 'w') as f:
        f.write(out_file)
        f.close()


# get scores
def tes_score(true, tfpn):
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for j in range(len(true)):
        if tfpn[j] == 0.0 and true[j] == 0.0:
            tp += 1
        if tfpn[j] == 1.0 and true[j] == 1.0:
            tn += 1
        if tfpn[j] == 0.0 and true[j] == 1.0:
            fp += 1
        if tfpn[j] == 1.0 and true[j] == 0.0:
            fn += 1
    acc = (tp + tn) / (tp + tn + fp + fn)
    if (tp + fn) == 0 or (tn + fp) == 0 or (tp + fp) == 0 or (tn + fn) == 0:
        return tp, tn, fp, fn, float('%.3f' % acc), float('%.3f' % 0), float('%.3f' % 0), float('%.3f' % 0), float(
            '%.3f' % 0)
    else:
        sn = tp / (tp + fn)
        sp = tn / (tn + fp)
        ppv = tp / (tp + fp)
        mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tn + fn) * (tp + fn) * (tn + fp))
        return tp, tn, fp, fn, float('%.3f' % acc), float('%.3f' % sn), float('%.3f' % sp), float(
            '%.3f' % ppv), float('%.3f' % mcc)


# 聚类
def tes_cluster(evaluate_score, evaluate_key):
    t_index = []
    s_index = []
    acc_score = {}
    for i in range(len(evaluate_score)):
        acc_score[evaluate_key[i]] = float('%.3f' % evaluate_score[i][4])
    for key in evaluate_key:
        t_tp = key.split('s')[0]
        s_tp = 's' + key.split('s')[1]
        if t_tp not in t_index:
            t_index.append(t_tp)
        if s_tp not in s_index:
            s_index.append(s_tp)
    cluster_t = visual_box(len(t_index))
    cluster_s = visual_box(len(s_index))
    t_number = visual_box(len(t_index))
    s_number = visual_box(len(s_index))
    for key in acc_score:
        t_tp = key.split('s')[0]
        s_tp = 's' + key.split('s')[1]
        cluster_t[t_index.index(t_tp)] += acc_score[key]
        cluster_s[s_index.index(s_tp)] += acc_score[key]
        t_number[t_index.index(t_tp)] += 1
        s_number[s_index.index(s_tp)] += 1
    for i in range(len(cluster_t)):
        cluster_t[i] = float('%.4f' % (cluster_t[i] / t_number[i]))
    for i in range(len(cluster_s)):
        cluster_s[i] = float('%.4f' % (cluster_s[i] / s_number[i]))
    return [cluster_t, cluster_s, t_index, s_index]


# document train
def tes_train_d(file, c_number, gamma, out, now_path):
    y, x = svm_train(file, float(c_number), float(gamma), os.path.join(now_path, out))
    visual_block(y)
    visual_block(x)


# folder train
def tes_train_f(folder, cg_path, out, now_path):
    file_path = os.path.join(now_path, folder)
    if out not in os.listdir(now_path):
        out_folder = os.path.join(now_path, out)
        os.makedirs(out_folder)
    else:
        out_folder = os.path.join(now_path, out)
    cg_box = read_hys(cg_path)
    start_e = 0
    for eachfile in os.listdir(file_path):
        start_e += 1
        easy_time(start_e, len(os.listdir(file_path)))
        out_path = os.path.join(out_folder, eachfile + '.model')
        y, x = svm_train(os.path.join(folder, eachfile), float(cg_box[eachfile][0]),
                         float(cg_box[eachfile][-1].strip('\n')), out_path)
        visual_block(y)
        visual_block(x)


# document evaluate
def tes_eval_d(file, c_number, gamma, crossv, out, now_path):
    print('')
    file = os.path.join(now_path, file)
    out_path = os.path.join(now_path, out + '-eval.txt')
    test_label, predict_label = svm_evaluate(file, float(c_number), float(gamma), int(crossv))
    score_line = tes_score(test_label, predict_label)
    line_num = ''
    for li in range(len(score_line)):
        each_num = score_line[li]
        if li < 4:
            line_num += str(int(each_num)) + '\t'
        else:
            line_num += str('%.3f' % each_num) + '\t'
    lines = 'Model Evaluation\n\ntp\ttn\tfp\tfn\tacc\tsn\tsp\tppv\tmcc\n' + line_num[:-1]
    lines += '\n\n' + read_intro('document_evaluate')
    with open(out_path, 'w') as f1:
        f1.write(lines)
        f1.close()


def tes_eval_f(folder, cg_path, crossv, out, now_path):
    print('')
    in_folder = os.path.join(now_path, folder)
    if out not in os.listdir(now_path):
        out_folder = os.path.join(now_path, out)
        os.makedirs(out_folder)
    cg_box = read_hys(cg_path)
    start_e = 0
    evaluate_score = []
    evaluate_key = []
    for eachfile in os.listdir(in_folder):
        start_e += 1
        easy_time(start_e, len(os.listdir(in_folder)))
        test_label, predict_label = svm_evaluate(os.path.join(folder, eachfile), float(cg_box[eachfile][0]),
                                                 float(cg_box[eachfile][-1].strip('\n')), int(crossv))
        score_line = tes_score(test_label, predict_label)
        evaluate_score.append(score_line)
        evaluate_key.append(eachfile.split('_')[0])
    # 输出文件
    out_lines = ''
    for i in range(len(evaluate_score)):
        line = evaluate_score[i]
        mid_line = evaluate_key[i]
        for j in line:
            mid_line += ',' + str(j)
        out_lines += mid_line + '\n'
    out_lines = 'index,tp,tn,fp,fn,acc,sn,sp,ppv,mcc\n' + out_lines
    with open(os.path.join(out, 'Features_eval.csv'), 'w') as f2:
        f2.write(out_lines)
    # 聚类
    cluster = tes_cluster(evaluate_score, evaluate_key)
    cluster_t = cluster[0]
    cluster_s = cluster[1]
    t_index = cluster[2]
    s_index = cluster[3]
    plot_eval_main(cluster_t, t_index, out, cluster_s, s_index, evaluate_score, evaluate_key)


def tes_search_d(file, now_path):
    best_c, best_g = svm_grid(os.path.join(now_path, file))
    print('\n>>>C_numbr: ' + str(best_c) + '\tGamma: ' + str(best_g) + '\n')


def tes_search_f(folder, now_path):
    print('')
    search_path = os.path.join(now_path, folder)
    key = ''
    start_e = 0
    for eachfile in os.listdir(search_path):
        start_e += 1
        easy_time(start_e, len(os.listdir(search_path)))
        best_c, best_g = svm_grid(os.path.join(search_path, eachfile))
        key += eachfile + '\tC_numbr: ' + str(best_c) + '\tGamma: ' + str(best_g) + '\n'
    with open(os.path.join(now_path, 'Hyperparameter.txt'), 'w') as f2:
        f2.write(key)


def tes_predict(file, model, out, now_path):
    test_label, predict_label = svm_predict(os.path.join(now_path, file), os.path.join(now_path, model), 5)
    result = 'True,Predict\n'
    for i in range(len(test_label)):
        result += str(test_label[i]) + ',' + str(predict_label[i]) + '\n'
    svm_acc(test_label, predict_label)
    with open(os.path.join(now_path, out) + '_predict.csv', 'w', encoding="utf-8") as f:
        f.write(result)
