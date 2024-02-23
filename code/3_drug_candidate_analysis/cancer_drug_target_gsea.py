##This code is available in python3##

import pandas as pd
import gseapy as gp
import os


PATH = "../../result/propagation"
out_PATH = "../../result/drug_result"
network = os.listdir(PATH)

print(network)
cancer_list = list(pd.read_csv("../../data/cancer_list.csv")["Type"])
cancer_list.remove("ESCA") ## no mapping DGA in ESCA

for net in network:
    par = pd.read_csv(PATH + net +"/pagerank_result.csv", index_col="Gene")
    par.index = par.index.astype(str)
    re = []
    for ty in cancer_list:
        print(ty)
        tpar = par[ty].sort_values(ascending=False)
        pre_res = gp.prerank(tpar, "../../data/drug/fda_cancer_drug.gmt",
                             outdir="../../result/gsea_supple/",
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
    re.to_csv(out_PATH + net +"/"+"fda_drug_target_gsea.csv")
        
        
