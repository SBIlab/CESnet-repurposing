##This code is available in python3##

import pandas as pd
import gseapy as gp
import os

PATH = "../../result/propagation"
network = os.listdir(PATH)
print(network)
GENESET = ["Reactome"]



for net in network[:12]:
    print(net)
    par = pd.read_csv(PATH + net +"/pagerank_result.csv", index_col="Gene")
    par.index = par.index.astype(str)
    par = par.dropna(how="all", axis=1)

    for gs in GENESET:

        try: os.mkdir(PATH + net + "/Reatome")
        except: pass
        re = []

        for ty in par.columns:
            tpar = par[ty].sort_values(ascending=False)
            pre_res = gp.prerank(tpar, "../../data/Reactome.v7.0.symbols.gmt",
                                 outdir="../../result/gsea_supple",
                                 permutation_num=1000, no_plot=True, processes=20, seed=10)
            terms = pre_res.res2d.index
            a = pd.DataFrame(pre_res.res2d)
            a = a.sort_values(by=["nes"], ascending=False)
            re.append(pd.DataFrame(a.iloc[0,:]).T)
            a.to_csv(PATH + net + "/Reactome/"+ty+".txt", sep="\t")

        
