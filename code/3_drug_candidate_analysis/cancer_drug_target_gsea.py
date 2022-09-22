##This code is available in python3##

import pandas as pd
import gseapy as gp
import os


df = pd.read_csv("/home/kwang/PROJECT/Cancer_Driver/data/pandrug/cancer_drug_all_target_new.csv", index_col="Type")
PATH = "../../result/propagation"
out_PATH = "../../result/drug_result"
network = os.listdir(PATH)

print(network)

for net in network:

    par = pd.read_csv(PATH + net +"/pagerank_result.csv", index_col="Gene")
    par.index = par.index.astype(str)
    re = []
    ttar = list(set(par.columns) & set(df.index))
    #ttar = ["BRCA"]
    for ty in ttar:
        print(ty)
        tpar = par[ty].sort_values(ascending=False)
        pre_res = gp.prerank(tpar, "/home/kwang/PROJECT/Cancer_Driver/data/pandrug/cancer_all_target_new.gmt",
                             outdir="/home/kwang/PROJECT/Cancer_Driver/result/pagerank/gsea/",
                             permutation_num=1000, no_plot=True, processes=3, min_size = 15, max_size = 1000, seed=10)
        terms = pre_res.res2d.index
        a = pd.DataFrame(pre_res.res2d)
        if ty in a.index:
            re.append(pd.DataFrame(pd.DataFrame(pre_res.res2d).loc[ty,:]).T)
    re = pd.concat(re)
    re.index.name = "Type"
    re = re.sort_index()
    try: os.mkdir(out_PATH + net)
    except: pass
    re.to_csv(out_PATH + net +"/"+"drug_target_gsea_d_para.csv")
        
        
