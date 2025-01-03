import pandas as pd
import numpy as np
import os
from scipy import stats

#ddf = pd.read_hdf("../../data/fitness_pandrugs_all_parsing_new.h5")
ddf = pd.read_table("../../data/drug/drug_info.txt", index_col="Standard Drug Name")
d_mat = {}
d_matx = {}



PATH = "../../result/propagation/"
out_PATH = "../../result/drug_result/"
network = os.listdir(PATH)
print network

for net in network:
    print net
    try: os.mkdir(out_PATH + net)
    except: pass
    #print dd
    prank = pd.read_csv(PATH + net +"/pagerank_result.csv", index_col="Gene")
    cancer_list = list(pd.read_csv("../../data/cancer_list.csv")["Type"])

    for ca in cancer_list:
        dval_list = []
        dmax_list = []
        #print ca
        for dg in ddf.index:
            tg = prank.loc[prank.index.isin(ddf.loc[dg,"Gene(s)"].split("|")), ca]
            d_val = np.sqrt(np.mean(tg**2))
            d_max = tg.max()
            dval_list.append(d_val)
        d_mat[ca] = dval_list
    d_mat = pd.DataFrame(d_mat, index = ddf.index)
    d_mat.index.name = "Drug"

    try: os.mkdir(out_PATH + net)
    except: pass
    d_mat.to_csv(out_PATH + net + "/drug_cancer_tc_score.csv")

    


