import pandas as pd
import numpy as np
import os
import random
from multiprocessing import Process, Queue
from scipy import stats
import pickle
from statsmodels.stats.multitest import multipletests

#ddf = pd.read_csv("../../data/drug/pandrug/fit/fitness_pandrugs_all_parsing_new_updated.csv", index_col="Standard Drug Name")
ddf = pd.read_table("../../data/drug/drug_info.txt", index_col="Standard Drug Name")


d_mat = {}
d_matx = {}

#df = pd.read_csv("../../data/drug/pandrug/fit/cancer_drug_all_target_new.csv", index_col="Type")
df = pd.read_csv("../../data/cancer_list.csv")
df.index = df["Type"]
PATH = "../../result/propagation/"
out_PATH = "../../result/drug_result/"
network = os.listdir(PATH)
#network = map(lambda x: x[:-19], network)
#network = network[4:]
print network

perm_num = 1000


def cal_rms(target):
    target = np.array(target)
    return np.sqrt(np.mean(target**2))




def perm_value(q_task, output):
    cancer_list = q_task.get()
    d_mat = {}
    for ca in cancer_list:
        dval_list = []
        for dg in ddf.index:
            target_num = len(ddf.loc[dg,"Gene(s)"].split("|"))
            tg = prank.loc[prank.index.isin(ddf.loc[dg,"Gene(s)"].split("|")), ca]
            if len(tg) > 0:
                rand_tg = []
                for rns in list(random_array[dg]):
                    rand_tg.append(cal_rms(prank.loc[rns, ca]))
                d_val = cal_rms(tg)
                r_mean = np.mean(rand_tg)
                r_stv = np.std(rand_tg)
                zscore = float(d_val - r_mean)/r_stv
                pval = stats.norm.sf(zscore)
                dval_list.append(pval)
            else:
                dval_list.append(np.nan)
        d_mat[ca] = dval_list
    d_mat = pd.DataFrame(d_mat, ddf.index)

    output.put(d_mat)

for net in network:
    print net
    try: os.mkdir(out_PATH + net)
    except: pass

    input_list  = []
    result = []
    
    prank = pd.read_csv(PATH + net +"/pagerank_result.csv", index_col="Gene")
    with open("../../data/drug/random_node/" + net + "_degree_ctrl.pickle","rb") as f:
        random_array = pickle.load(f)
    procs = []
    q_tasks = Queue()
    q_output = Queue()
    NUM_procs = 8

    if NUM_procs < len(df.index):
        for num in range(NUM_procs):
            input_list.append([])
        for num in range(len(df.index)):
            input_list[num % NUM_procs].append(df.index[num])
    else:
        NUM_procs = len(df.index)
        for num in range(NUM_procs):
            input_list.append([])
        for num in range(len(df.index)):
            input_list[num].append(df.index[num])

    for i in range(0, NUM_procs):
        procs.append(Process(target=perm_value, args=(q_tasks, q_output)))
    for p in procs:
        p.start()
    for j in range(len(procs)):
        q_tasks.put(input_list[j])
    for p in procs:
        result.append(q_output.get())
    for p in procs:
        p.join()
    d_mat = pd.concat(result, axis=1)
    d_mat.index.name = "Drug"
        #os.system("cp " + out_PATH + net +"/drug_cancer_pagerank_"+str(dd)+"_rms_significant_pvalue_all_new.csv "+ out_PATH + net +"/old/")
    d_mat.to_csv(out_PATH + net +"/drug_tc_score_pvalue.csv")
    fdr_mat = pd.DataFrame(np.full((len(d_mat.index), len(d_mat.columns)), np.nan), index=d_mat.index, columns=d_mat.columns)
    f_mat = {}
    for cc in d_mat.columns:
        f_mat[cc] = [d_mat[cc].dropna().index, multipletests(d_mat[cc].dropna(), alpha=0.05, method='fdr_bh')[1]]
    for cc in d_mat.columns:
        fdr_mat.loc[f_mat[cc][0],cc] = f_mat[cc][1]
    fdr_mat.index.name = "Drug"
    fdr_mat.to_csv(out_PATH + net +"/drug_tc_score_fdr.csv")

