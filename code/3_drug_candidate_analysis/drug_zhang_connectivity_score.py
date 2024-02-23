import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse
from multiprocessing import Pool
import os

def cal_zhang_score(pt):
    cmap = cdata[pt]
    cred = creds
    
    cmap_rank = np.abs(cmap).rank() * np.sign(cmap)
    cred_rank = np.abs(cred["fold"]).rank() * np.sign(cred["fold"])
    tar = list(set(cmap.index) & set(cred.index))
    crs = np.dot(cmap_rank.loc[tar], cred_rank.loc[tar])
    cmax = np.dot(np.array(range(len(cmap) - len(tar) + 1, len(cmap) + 1)),
                  np.array(range(len(cred) - len(tar) + 1, len(cred) + 1)))

    cscore = crs / float(cmax)

    return (pt, cscore)



gene_info = pd.read_table("../../data/cmap_phase2/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", dtype=str)
landmark_gene_row_ids = gene_info["pr_gene_id"][gene_info["pr_is_lm"] == "1"]
landmark_only_gctoo = parse("../../data/cmap_phase2/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx", rid = landmark_gene_row_ids)
creed = os.listdir("../../data/creeds/deg/")
gene_info = gene_info.set_index("pr_gene_id")
cdata = landmark_only_gctoo.data_df
cdata.index = gene_info.loc[cdata.index,"pr_gene_symbol"]


NUM_procs = 15
result = {}
for cr in ['dz:552']:
    creds = pd.read_table("../../data/creeds/deg/"+cr+".txt",names=["Gene", "fold"])
    creds = creds.set_index("Gene")
    print creds
    pool = Pool(processes=NUM_procs)
    re = pool.map(cal_zhang_score, landmark_only_gctoo.data_df.columns)
    re = pd.Series(dict(re))
    name = cr
    result[name] = re

result = pd.DataFrame(result)
result.index.name = "sig_id"
result.to_csv("../../result/cmap_score/cmap_zhang_score.csv")



