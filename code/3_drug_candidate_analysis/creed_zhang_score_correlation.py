import pandas as pd
import numpy as np
from scipy import stats
import os


PATH = "../../result/drug_result/"
creed_info = pd.read_csv("../../data/creeds/signature_meta.csv", index_col="id")
ddf = pd.read_table("../../data/drug/drug_info.txt")
ddf = ddf.dropna(subset=["ID"]).set_index("ID")
ddf = pd.DataFrame(ddf["Standard Drug Name"], index =ddf.index)
ddf.columns = ["drug"]
cmap = pd.read_table("../../data/cmap/cmap_pubchempy_map.txt")
cmap.columns = ["drug","ID"]
cmap = cmap.set_index("ID")
#dmap = pd.read_table("../../data/pandrug_id_map.txt")
#dmap.columns = ["drug","ID"]
dmap = ddf
#dmap = dmap.set_index("ID")
sinfo = pd.read_table("../../data/cmap/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt", index_col="sig_id")
cinfo = pd.read_table("../../data/cmap/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt", index_col="cell_id")
cinfo = cinfo.loc[cinfo["sample_type"]=="tumor",]
cinfo = cinfo.loc[cinfo["modification"]=="-666",]

df = pd.read_csv("../../result/cmap_score/cmap_zhang_score.csv", index_col="sig_id")
t2c = {"COADREAD":["large intestine"]}
cmap = cmap.loc[cmap.index.drop_duplicates(keep=False)]
dmap = dmap.loc[dmap.index.drop_duplicates(keep=False)]


NETWORK = os.listdir( "../../result/drug_result/")
NETWORK = ['co-essentiality']
target_drug = list(set(dmap.index) & set(cmap.index))
zs = {}
l1000 = {}
print NETWORK
can = "COADREAD"
cc = list(df.columns)[0]
for di in target_drug:
    max_list = []
    for cell in cinfo.loc[cinfo["primary_site"].isin(t2c[can]),].index:
        cell_sig = sinfo.loc[(sinfo["cell_id"] == cell) & (sinfo["pert_iname"] == cmap.loc[di, "drug"]),].index
        if len(cell_sig) > 0:
            cell_max = -np.min(df.loc[cell_sig, cc])
            max_list.append(cell_max)
    if len(max_list)>0:
        max_zhang = np.median(max_list)
        zs[di] = max_zhang
    else: pass
l1000[cc] = pd.Series(zs)




'''
cancer_type = df.columns
for cc in cancer_type:
    print cc
    zs = {}
    can = t2t[creed_info.loc[cc,"disease_name"]]
    if can == "COADREAD":
        for di in target_drug:
            max_list = []
            for cell in cinfo.loc[cinfo["primary_site"].isin(t2c[can]),].index:
                cell_sig = sinfo.loc[(sinfo["cell_id"] == cell) & (sinfo["pert_iname"] == cmap.loc[di, "drug"]),].index
                if len(cell_sig) > 0:
                    cell_max = -np.min(df.loc[cell_sig, cc])
                    max_list.append(cell_max)
            #if len(target_sig)>0:
            if len(max_list)>0:
                #max_zhang = -np.median(df.loc[target_sig, cc])
                max_zhang = np.median(max_list)
                zs[di] = max_zhang
            else: pass
        l1000[cc] = pd.Series(zs)
    else: pass
'''

re_s = {}
re_p = {}
TYPE = l1000.keys()
#TYPE = ["colorectal cancer_Colon"]
for net in NETWORK:
    print net
    ps = pd.read_csv(PATH + net +"/drug_cancer_tc_score.csv",
                     index_col="Drug")
    spearman = []
    pearson = []
    for cc in TYPE:
        zha_score = l1000[cc].dropna()
        net_score = ps.loc[dmap.loc[zha_score.index,"drug"], can]
        #net_score = ps.loc[dmap.loc[zha_score.index, "drug"], creed_info.loc[cc, "Type"]]
        net_score.index = zha_score.index
        net_score = net_score.dropna()
        zha_score = zha_score.loc[net_score.index,]
        spearman.append(stats.spearmanr(zha_score, net_score)[0])
    re_s[net] = spearman

re_s = pd.DataFrame(re_s, index=l1000.keys())
re_s.to_csv("../../result/drug_result/l1000_spearman_zhang_score_cell_median.csv")
