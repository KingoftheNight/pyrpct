# import packages
import os
from Read import read_intlen_hys, read_intro
from SVM import svm_train, svm_predict
import shutil
import math
from Visual import visual_block


# 排序
def intlen_sort(acc_value):
    arr = []
    for i in acc_value:
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


# 特征评估文件读取
def intlen_evaluate(eval_file, member, train_fs, predict_fs):
    with open(eval_file, 'r') as ef:
        data = ef.readlines()
        ef.close()
    eval_key = []
    acc_value = []
    for each in data[1:]:
        eachline = each.split(',')
        eval_key.append(eachline[0])
        acc_value.append(float(eachline[5]))
    acc_sort = intlen_sort(acc_value)
    t_fs_path = []
    p_fs_path = []
    fs_key = []
    for i in range(member):
        t_fs_path.append(os.path.join(train_fs, eval_key[acc_sort[i]] + '_rpct.fs'))
        p_fs_path.append(os.path.join(predict_fs, eval_key[acc_sort[i]] + '_rpct.fs'))
        fs_key.append(eval_key[acc_sort[i]] + '_rpct.fs')
    return [t_fs_path, p_fs_path, fs_key]


# 训练模型并复制预测特征文件
def intlen_train(t_fs_path, p_fs_path, fs_key, cg_box, now_path):
    if 'Intlen_predict' not in os.listdir(now_path):
        intlen_path = os.path.join(now_path, 'Intlen_predict')
        os.makedirs(intlen_path)
    else:
        intlen_path = os.path.join(now_path, 'Intlen_predict')
    test_model = []
    test_features = []
    for i in range(len(t_fs_path)):
        each_t_fs = t_fs_path[i]
        each_t_ms = os.path.join(intlen_path, fs_key[i] + '.model')
        test_model.append(each_t_ms)
        each_p_fs = p_fs_path[i]
        each_p_ms = os.path.join(intlen_path, fs_key[i])
        test_features.append(each_p_ms)
        each_c_num = cg_box[i][0]
        each_g_num = cg_box[i][1]
        y, x = svm_train(each_t_fs, each_c_num, each_g_num, each_t_ms)
        visual_block(y)
        visual_block(x)
        shutil.copyfile(each_p_fs, each_p_ms)
    return test_model, test_features


# 模型预测
def intlen_predict(test_model, test_features):
    predict_value = []
    y_p = []
    for i in range(len(test_model)):
        each_model = test_model[i]
        each_feature = test_features[i]
        y_p, p_label = svm_predict(each_feature, each_model, 5)
        predict_value.append(p_label)
    true_value = y_p
    return true_value, predict_value


# 多数投票
def intlen_vote(true_value, predict_value):
    out_value = []
    for i in range(len(true_value)):
        mid_value = []
        for each_value in predict_value:
            mid_value.append(each_value[i])
        number_0 = 0
        number_1 = 0
        for j in mid_value:
            if j == 0.0:
                number_0 += 1
            else:
                number_1 += 1
        if number_0 > number_1:
            out_value.append(0.0)
        else:
            out_value.append(1.0)
    # 集合模型评估
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for j in range(len(true_value)):
        if true_value[j] == 0.0 and out_value[j] == 0.0:
            tp += 1
        if true_value[j] == 1.0 and out_value[j] == 1.0:
            tn += 1
        if true_value[j] == 0.0 and out_value[j] == 1.0:
            fp += 1
        if true_value[j] == 1.0 and out_value[j] == 0.0:
            fn += 1
    acc = (tp + tn) / (tp + tn + fp + fn)
    sn = tp / (tp + fn)
    sp = tn / (tn + fp)
    ppv = tp / (tp + fp)
    mcc = (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tn + fn) * (tp + fn) * (tn + fp))
    # 结果输出
    out_line = str(tp) + '\t' + str(tn) + '\t' + str(fp) + '\t' + str(fn) + '\t' + str('%.3f' % acc) + '\t' + str(
        '%.3f' % sn) + '\t' + str('%.3f' % sp) + '\t' + str('%.3f' % ppv) + '\t' + str('%.3f' % mcc)
    lines = 'Model Evaluation\n\ntp\ttn\tfp\tfn\tacc\tsn\tsp\tppv\tmcc\n' + out_line
    lines += '\n\n' + read_intro('document_evaluate')
    return lines


# 模型集合投票主程序
def intlen_main(train_fs, predict_fs, eval_file, cg_file, member, now_path):
    # 参数格式修改
    train_fs = os.path.join(now_path, train_fs)
    predict_fs = os.path.join(now_path, predict_fs)
    eval_file = os.path.join(now_path, eval_file)
    cg_file = os.path.join(now_path, cg_file)
    member = int(member)
    # 读取评估文件
    features_files_path = intlen_evaluate(eval_file, member, train_fs, predict_fs)
    t_fs_path = features_files_path[0]
    p_fs_path = features_files_path[1]
    fs_key = features_files_path[2]
    # 读取超参数文件
    cg_box = read_intlen_hys(fs_key, cg_file)
    # 训练模型
    test_model, test_features = intlen_train(t_fs_path, p_fs_path, fs_key, cg_box, now_path)
    # 模型预测
    true_value, predict_value = intlen_predict(test_model, test_features)
    # 多数投票
    out_content = intlen_vote(true_value, predict_value)
    with open(os.path.join(now_path, 'Integrated_learning_' + str(member) + '.txt'), 'w') as lf:
        lf.write(out_content)
