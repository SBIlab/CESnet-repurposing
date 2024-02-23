import pandas as pd
import numpy as np
from scipy import stats
import os
import pickle


ddf = pd.read_table("../../data/drug/drug_info.txt", index_col="Standard Drug Name")
ddf = ddf.dropna(subset=['ID'])
prism_map = pd.read_table("../../data/drug/prism_id_map.txt", index_col="drug_name")
#pandrug_map = pd.read_table("../../data/pandrug_id_map.txt", index_col="drug_name")
pandrug_map = pd.DataFrame(ddf["ID"])
rpo = pd.read_hdf("../../data/drug/prism_ic50_matrix.h5")
rpo = rpo[prism_map.index]
rpo = rpo[rpo.columns.drop_duplicates(keep=False)]
rpo.columns = prism_map["ID"]
rpo = rpo.replace([np.inf, -np.inf], np.nan)
rpo = -np.log2(rpo)
rpo = rpo.dropna(axis="columns", thresh=48)
rpo = rpo.replace([np.inf], 40)


with open("../../data/BR20q2_celltype_mapping.pickle", "rb") as f:
    ctype = pickle.load(f)


d_para = [0.85]
df = pd.read_csv("../../data/drug/fda_cancer_drug.csv", index_col="Type")


PATH = "../../data/network/"
out_PATH = "../../result/drug_result/"

network = os.listdir(out_PATH)
network = filter(lambda x: not x=="old", network)

print network
for net in network:
    print net
    try: os.mkdir(out_PATH + net) # os.mkdir(out_PATH + net +"/repurposing/new/")
    except: pass
    d_mat = pd.read_csv(out_PATH + net +"/drug_cancer_tc_score.csv",index_col="Drug")
    d_mat = d_mat.loc[pandrug_map.index,]
    d_mat.index = pandrug_map["ID"]
    d_mat = d_mat.loc[d_mat.index.drop_duplicates(keep=False)]
    out_index = []
    result = []
    SP = []
    SPP = []
    dmc = list(set(d_mat.columns)&set(ctype.keys()))
    #dmc.remove("CESC")

    for i, ca in enumerate(dmc):
        approve =[]
        if ca in df.index:
            tmp_target = df.loc[ca, "drug"].split("|")
            for tt in tmp_target:
                if tt in pandrug_map.index:
                    approve.append(pandrug_map.loc[tt, "ID"])
            approve = list(set(approve) & set(d_mat.index))
        other = d_mat.loc[~d_mat.index.isin(approve), ca].dropna()
        other = other.loc[other.index.isin(rpo.columns)]
        other = other.index
        other_value = rpo.loc[rpo.index.isin(ctype[ca]), rpo.columns.isin(other)].median()
        other_value = other_value.dropna()

        #out_index.append(ca)
        Y = list(other_value)
        X = list(d_mat.loc[other_value.index, ca])
        if len(X)>0:
            out_index.append(ca)
            SP.append(stats.spearmanr(X,Y)[0])
            SPP.append(stats.spearmanr(X,Y)[1])
    #output = pd.DataFrame([SP, SPP], columns=["spearmanr", "spearman-pval"], index=out_index)
    output = pd.DataFrame({"spearmanr":SP, "spearman-pval":SPP} ,index=out_index)
    output.index.name = "Type"
    output = output.sort_index(axis=0)
    output.to_csv(out_PATH + net + "/drug_tc_score_ic50_spearmanr_pubchem_median.csv")
