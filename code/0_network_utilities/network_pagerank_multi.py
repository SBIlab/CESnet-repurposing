import pandas as pd
import numpy as np
import networkx as nx
from multiprocessing import Process, Queue
import os

def read_module(path, target):
    mod = {}
    f = open(path, "r")
    for line in f.xreadlines():
        line = line.strip().split("\t")
        if line[0] in target:
            mod[line[0]] = line[2:]
    return mod



def page_rank(q_task, output):
    G, mod, dp = q_task.get()
    result = []
    node_list=list(G.nodes())
    for k in mod.keys():
        p_vec = {k: v for k, v in dict(zip(mod[k], [1]*len(mod[k]))).items() if k in node_list}
        if len(p_vec) > 0:
            re = nx.pagerank_scipy(G, personalization=p_vec, alpha=dp)
            result.append(pd.DataFrame.from_dict(re, orient="index", columns=[k]).sort_index())
        else:
            result.append(pd.DataFrame({k: [np.nan]*len(node_list)}, index=node_list).sort_index())
    result = pd.concat(result, axis=1)
    output.put(result)


target_type = []
with open("../../data/target_cancer_type.txt") as f:
    for line in f.xreadlines():
        target_type.append(line.strip())


dp = 0.85
PATH = "../../data/network/"
network = os.listdir(PATH)
file_list = network
print file_list


for fk in range(len(file_list)):
    print fk
    net = file_list[len(file_list)-fk-1]
    print net[:-4]
    with open(PATH + net, "rb") as f:
        f.readline()
        G = nx.read_edgelist(f, delimiter="\t", nodetype=str, data=(('weight',float),))

    # read driver gene
    modul = read_module("../../data/driver/Bailey_driver_symbol.gmt", target_type)

    #multiprocessing
    print "Now start Multi-Processing!"
    output = []
    procs = []
    q_tasks = Queue()
    q_output = Queue()
    result = []
    NUM_procs = 5
    mk = modul.keys()
    input_list = []

    if len(mk) > NUM_procs:
        for num in range(NUM_procs):
            input_list.append(dict())
        for n in range(len(mk)):
            input_list[n%NUM_procs][mk[n]] = modul[mk[n]]
    else:
        NUM_procs = len(mk)
        for num in range(NUM_procs):
            input_list.append(dict())
        for n in range(len(mk)):
            input_list[n%NUM_procs][mk[n]] = modul[mk[n]]

    for i in range(0, NUM_procs):
        procs.append(Process(target=page_rank, args=(q_tasks, q_output)))
    for p in procs:
        p.start()
    for j in range(len(procs)):
        q_tasks.put([G, input_list[j], dp])
    for p in procs:
        result.append(q_output.get())
    for p in procs:
        p.join()
    outpath = "../../result/propagation/" + file_list[fk][:-4]
    try: os.mkdir(outpath)
    except:pass
    result = pd.concat(result, axis=1)
    result.index.name = "Gene"
    result.to_csv(outpath + "/pagerank_result.csv")
    del G




