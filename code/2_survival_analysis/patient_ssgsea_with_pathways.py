##This code is available in python3##

import pandas as pd
import os
import gseapy as gp

PATH = "../../data/TCGA_patient/"
GENESET = ["Reactome"]
for gs in GENESET:
    TYPE = os.listdir(PATH)
    TYPE.remove("Reactome")
    print(TYPE)
    try: os.mkdir("../../data/" + gs)
    except: pass
    for ty in TYPE:
        print(ty[5:])
        em = pd.read_table(PATH+ty+"/GeTMM_log2_rna_seq.txt", index_col="genes")
        ssres = gp.ssgsea(em, "../../data/Reactome.v7.0.ensembl.gmt",
                             outdir="../../result/pagerank/gsea/", sample_norm_method='rank',
                             permutation_num=0, processes=5,no_plot=True, min_size=5, seed=10)
        re = ssres.res2d.sort_index()
        print(re)
        re = pd.DataFrame(re)
        re.to_csv("../../data/Reactome/"+ty[5:] +".csv")

