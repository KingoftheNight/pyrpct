# import packages
import os
import time
import math
import copy
from Visual import visual_aa
res_path = os.path.dirname(__file__)


# 欧氏距离
def res_ojld(a, b):
    sq = 0
    for i in range(len(a)):
        sq += (a[i] - b[i]) * (a[i] - b[i])
    distance = math.sqrt(sq)
    return distance


# pre_c
def pre_r(euclidean_box, index):
    mid_list = []
    for key in euclidean_box:
        mid_list.append(euclidean_box[key])
    mid_list.sort()
    box = []
    for each_num in mid_list:
        for key in euclidean_box:
            if euclidean_box[key] == each_num:
                if key not in box:
                    box.append(key)
    pre_raac = []
    for m in mid_list:
        for key in euclidean_box:
            if euclidean_box[key] == m:
                if [key.split("&")[0], key.split("&")[1]] not in pre_raac:
                    pre_raac.append([key.split("&")[0], key.split("&")[1]])
    reduce_list = []
    aa_raac = copy.deepcopy(pre_raac)
    aa20 = copy.deepcopy(index)
    for i in aa_raac[:190]:
        if i[0] in aa20 and i[1] in aa20:
            aa20.remove(str(i[0]))
            aa20.remove(str(i[1]))
            aa20.append(i)
        else:
            p = 0
            q = 0
            if i[0] in aa20:
                aa20.remove(str(i[0]))
            if i[1] in aa20:
                aa20.remove(str(i[1]))
            for j in range(len(aa20)):
                if len(aa20[j]) == 1:
                    pass
                else:
                    if i[0] in aa20[j] or i[1] in aa20[j]:
                        p += 1
                        if p == 1:
                            if i[0] not in aa20[j]:
                                aa20[j].append(str(i[0]))
                            if i[1] not in aa20[j]:
                                aa20[j].append(str(i[1]))
                            q = copy.deepcopy(j)
                        else:
                            for k in aa20[j]:
                                if k not in aa20[q]:
                                    aa20[q].append(k)
                            aa20.remove(aa20[j])
                            break
        result = ""
        for amp in aa20:
            if len(amp) != 1:
                for poi in amp:
                    result += poi
                result += "-"
            else:
                result += amp + "-"
        result = result[:-1]
        if result not in reduce_list:
            reduce_list.append(str(result))
    return reduce_list


# 约化氨基酸
def res_reduce(data):
    index = visual_aa()
    name = data[0]
    data = data[1:]
    euclidean_box = {}
    print('>>>>>>>>>>>>>>>>>')
    m = 0
    for i in range(20):
        m += 1
        for j in range(m, 20):
            tube_a = []
            tube_b = []
            tube_a.append(float(data[i]))
            tube_b.append(float(data[j]))
            # 欧氏距离
            distance = res_ojld(tube_a, tube_b)
            print('\r>>>' + index[i] + "&" + index[j] + ' 欧式距离：' + str('%.4f' % distance), end='', flush=True)
            euclidean_box[index[i] + "&" + index[j]] = str('%.4f' % distance)
            time.sleep(0.01)
    print('>>>>>>>>>>>>>>>>>')
    final = []
    reduce_list = pre_r(euclidean_box, index)
    for y in reduce_list:
        result = name + " 1 size " + str(len(y.split("-"))) + " " + y
        final.insert(0, result)
    final = final[1:]
    return final


# res main
def res_main(resp_id):
    aaindex_path = os.path.join(os.path.join(res_path, 'aaindexDB'), 'AAindex.txt')
    with open(aaindex_path, 'r', encoding='UTF-8') as f:
        resp_data = f.readlines()
        f.close()
    request = ''
    for i in resp_data:
        if resp_id in i:
            request = i
    if len(request) == 0:
        print('We do not find the ID:' + resp_id)
        return
    print("\n约化中...")
    final = res_reduce(request.split('\t'))
    print("约化完成.")
    time.sleep(0.5)
    # 打印约化列表
    print('>>>>>>>>>>>>>>>>>')
    for j in final:
        time.sleep(0.05)
        print(j)
    print('>>>>>>>>>>>>>>>>>')
    print("结果打印完毕.")
