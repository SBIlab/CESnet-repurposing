import numpy as np
import pandas as pd
from scipy import stats
from sklearn import metrics
from multiprocessing import Process, Queue
import os
import networkx as nx

def cal_AUC(X, Y):
    fpr, tpr, thresholds = metrics.roc_curve(Y, X, pos_label=1)
    AUC = metrics.auc(fpr, tpr)
    # f1 = plt.plot(fpr, tpr)
    # f1.savefig()

    return AUC

def read_module(path):
    mod = {}
    f = open(path, "r")
    for line in f.xreadlines():
        line = line.strip().split("\t")
        mod[line[0]] = line[2:]
    return mod


def filtering(mod, g_list):
    for mm in mod.keys():
        check = filter(lambda x: x in g_list, mod[mm])
        mod[mm] = check

    return 0


###
def performance_check(q_task, output):
    m_dict, cor = q_task.get()
    mX = {}
    mY = {}
    for mk in m_dict.keys():
        md = m_dict[mk]
        Y = []  # (len(md) * [1]) + (len(rot.index)-len(md)) *[0]
        X = []
        for gg in cor.columns:
            if gg in md:
                Y.append(1)
                ta = md[:]
                ta.remove(gg)
                zscore = cor.loc[ta, gg].sum()
                bot = cor.loc[:,gg].sum()
                if bot > 0: X.append(zscore/bot)
                else: X.append(0.)

            else:
                Y.append(0)
                zscore = cor.loc[md, gg].sum()
                bot = cor.loc[:,gg].sum()
                if bot > 0: X.append(zscore/bot)
                else: X.append(0.)

        mX[mk] = X
        mY[mk] = Y
    mX = pd.DataFrame(mX, index=cor.columns)
    mY = pd.DataFrame(mY, index=cor.columns)
    mX.index.name = "Gene"
    mY.index.name = "Gene"

    output.put((mX, mY))



####read_module

module_list = ["Driver_Balley"]
PATH = "../../data/network/edge_last/"
out_PATH = "../../result/driver_prediction_guilt_by_assoc"
NETWORK = os.listdir(PATH)

print NETWORK


for net in NETWORK:
    print net
    with open(PATH + net, "rb") as f:
        f.readline()
        G = nx.read_edgelist(f, delimiter="\t", nodetype=str, data=(('weight',float),))
    adj = nx.adjacency_matrix(G, weight='weight')
    cor = pd.DataFrame(adj.toarray(), columns=list(G.nodes()), index=list(G.nodes()))
    roc_auc = []
    p_list = []
    ind = []
    # zmat = {}
    modul = read_module("../../data/driver/Bailey_driver_symbol.gmt")
    filtering(modul, cor.index)
    ##================================================================
    print('Now performing multiprocessing test..')

    input_list = []
    NUM_procs = 12
    nu = len(cor.columns) / NUM_procs
    for i in range(NUM_procs):
        if i != NUM_procs - 1:
            input_list.append(cor.iloc[:, nu * i:nu * (i + 1)])
        else:
            input_list.append(cor.iloc[:, nu * i:])

    output = []
    procs = []
    au = []
    pv = []
    result = []
    # z_re = []
    mX = []
    mY = []
    q_tasks = Queue()
    q_output = Queue()

    for i in range(0, NUM_procs):
        procs.append(Process(target=performance_check, args=(q_tasks, q_output)))

    for p in procs:
        p.start()
    for j in range(len(procs)):
        q_tasks.put([modul, input_list[j]])
    for p in procs:
        a, b = q_output.get()
        mX.append(a)
        mY.append(b)
    for p in procs:
        p.join()
    to_mX = pd.concat(mX)
    to_mY = pd.concat(mY)
    to_mX = to_mX.sort_index()
    to_mY = to_mY.sort_index()


    for cc in to_mX.columns:
        au.append(cal_AUC(to_mX[cc], to_mY[cc]))
    for cc in to_mX.columns:
        ann = []
        oth = []
        for i, a in enumerate(to_mY[cc]):
            if a > 0:
                ann.append(to_mX[cc][i])
            else:
                oth.append(to_mX[cc][i])
        try: pv.append(stats.mannwhitneyu(ann, oth, alternative="greater")[1])
        except: pv.append(np.nan)
    to_re = pd.DataFrame({"AUC": au, "p-value": pv}, index=to_mX.columns)
    siz = []
    for mi in to_re.index:
        siz.append(len(modul[mi]))
    to_re["size"] = siz
    to_re.index.name = "Type"
    to_re.sort_index()
    npath = out_PATH + net[:-4] +"/"
    try:
        os.mkdir(npath)
    except:
        pass
    #ROAUC result
    to_re.to_csv(npath + net[:-4] +"_modularity_ROAUC.csv")
    #X value result
    to_mX.to_csv(npath + net[:-4] +"_X_value.csv")
    #Y value result
    to_mY.to_csv(npath + net[:-4] +"_Y_value.csv")


# ==============================================
