import pandas as pd
import numpy as np
import os
from khpylib import performance
from multiprocessing import Process, Queue
import pickle

def cal_performance(q_task, output):
    net, cancer = q_task.get()
    t_result = {}
    for ca in cancer:
        roauc = []
        prauc = []
        t_roauc = []
        t_prauc = []
        for i in range(100):
            PATH = "../../result/ukin/cgc_cgc/a0.85r0.85/%s/%s/" % (net, ca)
            rPATH = '../../data/ukin/random_set/cgc_cgc/'
            ukin = pd.read_table(PATH + "guiderw_result_%s.txt" % i, names=["gene", "grw"])
            ukin = ukin.set_index("gene")
            with open(rPATH + 'random_set_%s.pickle' % i, 'rb') as f:
                rset = pickle.load(f)
            pos = rset["answer"]
            exc = list(set(cgc) - set(pos))
            tar = list(set(ukin.index) - set(exc))
            ukin = ukin.loc[tar,]
            ukin["rank"] = ukin["grw"].rank(method='average', ascending=False)

            label = []
            for gg in ukin.index:
                if gg in pos:
                    label.append(1)
                else:
                    label.append(0)
            ukin["labels"] = label
            t100 = ukin.loc[ukin["rank"] <= 100,]
            labels = ukin["labels"]
            predicts = ukin["grw"]
            roauc.append(performance.cal_ROAUC(predicts, labels))
            prauc.append(performance.cal_PRAUC(predicts, labels))
            t_roauc.append(performance.cal_ROAUC(t100["grw"], t100["labels"]))
            t_prauc.append(performance.cal_PRAUC(t100["grw"], t100["labels"]))
            temp = pd.DataFrame({'roauc':roauc, 'prauc':prauc, 't_roauc':t_roauc, 't_prauc':t_prauc})
            temp.index.name = 'number'
            temp.to_csv(PATH+"summary_result.csv")
        t_result[ca] = [roauc, prauc, t_roauc, t_prauc]
    output.put(t_result)


NETWORK = os.listdir("../../result/ukin/cgc_cgc/a0.85r0.85")
CANCER = os.listdir("../../result/ukin/cgc_cgc/a0.85r0.85/fit_knn_2.0_cor_clr_abs")
CANCER = filter(lambda x: not x.startswith('prior'), CANCER)
#CANCER.remove('KICH')
#CANCER.remove('PRAD')
with open("../../data/ukin/driver.txt",'rb') as f:
    cgc = map(lambda x: x.strip(), f.readlines())
roauc_=dict()
prauc_=dict()
t_roauc_=dict()
t_prauc_=dict()
print NETWORK
NUM_procs = len(CANCER)
NUM_procs = 15

input_list = []
for i in range(NUM_procs):
    input_list.append([])
for i in range(len(CANCER)):
    input_list[i%NUM_procs].append(CANCER[i])



for net in NETWORK:
    print net
    c_roauc = {}
    c_prauc = {}
    c_t_roauc = {}
    c_t_prauc = {}
    q_tasks = Queue()
    q_output = Queue()
    procs = []
    for i in range(0, NUM_procs):
        procs.append(Process(target=cal_performance, args=(q_tasks, q_output)))
    for p in procs:
        p.start()
    for j in range(len(procs)):
        q_tasks.put([net,input_list[j]])
    for p in procs:
        t_dict= q_output.get()
        for ca in t_dict.keys():
            roauc, prauc, t_roauc, t_prauc = t_dict[ca]
            c_t_roauc[ca] = np.mean(filter(lambda x: not np.isnan(x),t_roauc))
    for p in procs:
        p.join()
    t_roauc_[net] = pd.Series(c_t_roauc)
t_roauc_ = pd.DataFrame(t_roauc_, index=CANCER).to_csv("../../result/ukin/cgc_cgc_a0.85r0.85_top100_roauc.csv")











