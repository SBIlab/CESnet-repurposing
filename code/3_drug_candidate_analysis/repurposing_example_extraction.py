import pandas as pd
import numpy as np
from scipy import stats
import os
import pickle

rpo = pd.read_hdf("../../data/prism_ic50_matrix.h5")

with open("../../data/drug/prism/BR20q2_celltype_mapping.pickle", "rb") as f:
    ctype = pickle.load(f)

dr = {}
with open("../../data/driver/Bailey_driver_symbol.gmt", "rb") as f:
    for line in f.xreadlines():
        line = line.strip().split("\t")
        dr[line[0]] = line[2:]

ddf = pd.read_hdf("../../data/fitness_pandrugs_all_parsing_new.h5")

d_para = [0.85]
df = pd.read_csv("../../data/cancer_drug_all_target_new.csv", index_col="Type")
out_PATH = "../../result/drug_result/"
network = os.listdir("../../data/network")
TYPE = list(pd.read_csv("../../data/cancer_list.csv")["Type"])
print ctype.keys()

print network
print TYPE
for da in d_para:
    for i, ca in enumerate(TYPE):
        print ca
        result = {}
        for net in network:

            d_mat = pd.read_csv(out_PATH + net +"/drug_cancer_tc_score.csv",index_col="Drug")

            out = pd.Series(np.zeros(len(d_mat)), index = d_mat.index)
            target = list(set(df.loc[ca, "drug"].split("|")) & set(d_mat.index))
            other = d_mat.loc[~d_mat.index.isin(target), ca].dropna()
            cand = pd.Series(map(lambda x: 100.-stats.percentileofscore(other, x), other), index=other.index)

            out.loc[target] = -1
            out.loc[cand.index] += cand
            result[net] = out
        result = pd.DataFrame(result, index=d_mat.index)
        result["Source Drug Name"] = ddf.loc[result.index,"Source Drug Name"]
        result["Status"] = ddf.loc[result.index,"Status"]
        result["Pathology"] = ddf.loc[result.index,"Pathology"]
        result["Cancer(s)"] = ddf.loc[result.index,"Cancer(s)"]
        result.loc[result["Cancer(s)"]=="CANCER", "Pathology"] = "CANCER"
        del result["Cancer(s)"]
        check = []
        for rr in result.index:
            if result.loc[rr,"fit_knn_2.0_cor_clr_abs"] <= 10: ##9090 * 10%
                if len(filter(lambda x: x <= 10, result.loc[rr,network])) == 1:
                    check.append(1)
                else:
                    check.append(0)
            else:
                check.append(0)
        result["min_check"] = check
        result["target"] = ddf.loc[result.index,"Gene(s)"]
        non_driver_target = []
        for rr in result.index:
            ttp = set(ddf.loc[rr,"Gene(s)"].split("|"))
            ttp = list(ttp - set(dr[ca]))
            non_driver_target.append("|".join(ttp))
        result["ND_target"] = non_driver_target
        result.to_csv("../../result/repurposing/"+ca+".txt", sep="\t")


