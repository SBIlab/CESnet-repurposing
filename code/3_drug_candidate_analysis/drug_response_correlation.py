import pandas as pd
import numpy as np
from scipy import stats
import os
import pickle

prism_map = pd.read_table("../../data/prism_id_map.txt", index_col="drug_name")
pandrug_map = pd.read_table("../../data/pandrug_id_map.txt", index_col="drug_name")
rpo = pd.read_hdf("../../data/prism_ic50_matrix.h5")
rpo = rpo[prism_map.index]
rpo = rpo[rpo.columns.drop_duplicates(keep=False)]
rpo.columns = prism_map["ID"]
rpo = rpo.replace([np.inf, -np.inf], np.nan)
rpo = -np.log2(rpo)
rpo = rpo.dropna(axis="columns", thresh=48)
rpo = rpo.replace([np.inf], 40)


with open("../../data/BR20q2_celltype_mapping.pickle", "rb") as f:
    ctype = pickle.load(f)

ddf = pd.read_hdf("../../data/fitness_pandrugs_all_parsing_new.h5")

d_para = [0.85]
df = pd.read_csv("../../data/cancer_drug_all_target_new.csv", index_col="Type")


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
    par = pd.read_csv(PATH + net + "/pagerank_result.csv", index_col="Gene")
    out_index = []
    result = []
    SP = []
    SPP = []
    dmc = list(set(d_mat.columns)&set(ctype.keys()))
    dmc.remove("CESC")

    for i, ca in enumerate(dmc):
        tmp_target = df.loc[ca, "drug"].split("|")
        target = []
        for tt in tmp_target:
            if tt in pandrug_map.index:
                target.append(pandrug_map.loc[tt, "ID"])
        target = list(set(target) & set(d_mat.index))

        other = d_mat.loc[~d_mat.index.isin(target), ca].dropna()
        other = other.loc[other.index.isin(rpo.columns)]
        appo = target
        other = other.index
        ohter_value = rpo.loc[rpo.index.isin(ctype[ca]), rpo.columns.isin(other)].median()
        other_value = other_value.dropna()

        out_index.append(ca)
        Y = list(other_value)
        X = list(d_mat.loc[other_value.index, ca])
        if len(X)>0:
            tdf = pd.DataFrame({"X":np.log10(X), "Y": Y})
            SP.append(stats.spearmanr(X,Y)[0])
            SPP.append(stats.spearmanr(X,Y)[1])
    output = pd.DataFrame([SP, SPP], columns=["spearmanr", "spearman-pval"], index=out_index)
    output.index.name = "Type"
    output = output.sort_index(axis=0)
    output.to_csv(out_PATH + net + "/drug_tc_score_ic50_spearmanr_pubchem_median.csv")
