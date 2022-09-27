import networkx as nx
import pandas as pd
import numpy as np
import os, sys
# import random
import pickle

def page_rank(G, d_factor, p_vec):
    nodes= list(G.nodes())
    p_vec = {k: v for k, v in p_vec.items() if k in nodes}
    if len(p_vec) > 0:
        p_result = nx.pagerank_scipy(G, alpha=d_factor, personalization=p_vec)
    else:
        p_result = dict(zip(nodes,[np.nan]*len(nodes)))
    return pd.Series(p_result).sort_index()



alpha, gamma, network, cancer, num = sys.argv[1:]
#fit.txt, 0.5, 1.0, BRCA, 4
alpha = float(alpha)
gamma = float(gamma)
num = int(num)

try:
    os.mkdir("../../result/ukin/cgc_cgc/a%sr%s" %(alpha, gamma))
except:
    pass
try:
    os.mkdir("../../result/ukin/cgc_cgc/a%sr%s/%s" %(alpha, gamma, network[:-4]))
except:
    pass
try:
    os.mkdir("../../result/ukin/cgc_cgc/a%sr%s/%s/%s" %(alpha, gamma, network[:-4],cancer[:-4]))
except:
    pass


#################
#prior knowledge
#################

with open("../../data/network/" + network, "rb") as f:
    f.readline()
    G = nx.read_edgelist(f, delimiter="\t", nodetype=str, data=(('weight', float),))


#with open("../../data/ukin/driver.txt", "rb") as f:
#    cgc = map(lambda x: x.strip(),f.readlines())


if os.path.isfile("../../result/ukin/cgc_cgc/a%sr%s/%s/prior_result_%s.pickle" % (alpha, gamma, network[:-4],num)):
    with open("../../result/ukin/cgc_cgc/a%sr%s/%s/prior_result_%s.pickle"% (alpha, gamma, network[:-4],num)) as f:
        prior_result = pickle.load(f)


else:
    with open("../../data/ukin/random_set/cgc_cgc/random_set_%s.pickle" % num, "rb") as f:
        random_input = pickle.load(f)
    answer = random_input["answer"]
    prior = random_input["prior"]
    prior_result = page_rank(G, gamma, dict(zip(prior, [1] * len(prior))))
    with open("../../result/ukin/cgc_cgc/a%sr%s/%s/prior_result_%s.pickle" % (alpha, gamma, network[:-4],num),"wb") as fw:
        pickle.dump(prior_result, fw, protocol=pickle.HIGHEST_PROTOCOL)

'''
cgc=set(cgc) & set(G.nodes)
hidden = random.sample(cgc, 400)
left = cgc - set(hidden)
prior = list(random.sample(left, 20))
answer = list(hidden)
'''


if prior_result.isnull().any():
    guiderw_result = pd.Series([np.nan]*len(prior_result), index=prior_result.index).sort_index()
    guiderw_result.to_csv("../../result/ukin/cgc_cgc/a%sr%s/%s/%s/guiderw_result_%s.txt" % (alpha, gamma, network[:-4], cancer[:-4], num),sep="\t")

else:
    ###########################
    #new information
    ###########################

    # make new network
    D = nx.DiGraph()
    for nd in G.nodes():
        neibor = list(G.neighbors(nd))
        w_sum = prior_result.loc[neibor].sum()
        for ne in neibor:
            D.add_edge(nd, ne, weight=prior_result.loc[ne]/w_sum)

    d_node = D.nodes()
    #read_new_info
    new_info = {}
    with open("../../data/ukin/normalized_mutation/%s" % (cancer),"rb") as f:
        for line in f.xreadlines():
            line = line.strip().split("\t")
            if line[0] in d_node():
                new_info[line[0]] = float(line[1])

    #guided random walk
    guiderw_result = page_rank(D, alpha, new_info)

    guiderw_result.to_csv("../../result/ukin/cgc_cgc/a%sr%s/%s/%s/guiderw_result_%s.txt" %(alpha, gamma, network[:-4],cancer[:-4],num), sep="\t")



'''
with open("../../result/ukin/cgc_cgc/a%sr%s/%s/%s/answer_result_%s.txt" %(alpha, gamma, network[:-4],cancer[:-4],num),"wb") as fw:
    for aa in answer:
        fw.write(aa+"\n")
with open("../../result/ukin/cgc_cgc/a%sr%s/%s/%s/prior_result_%s.txt" %(alpha, gamma, network[:-4],cancer[:-4],num),"wb") as fw:
    for pr in prior:
        fw.write(pr+"\n")
'''






