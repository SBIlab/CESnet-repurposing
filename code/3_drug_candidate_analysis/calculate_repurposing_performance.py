import pandas as pd
import numpy as np
from sklearn import metrics
import os
from collections import defaultdict
from khpylib import performance

def get_TPFNTNFP_from_strip(actual_list, predicted_list):
    TP, FN, TN, FP = 0, 0, 0, 0
    for i in range(len(actual_list)):
        if actual_list[i] == 1:
            if predicted_list[i] == 1: TP += 1
            elif predicted_list[i] == 0: FN += 1
        elif actual_list[i] == 0:
            if predicted_list[i] == 1: FP += 1
            elif predicted_list[i] == 0: TN += 1
    return [TP, FN, TN, FP]

def F1_score(TPFNTNFP, cri=['num', 'list'][0]):
    TP, FN, TN, FP = TPFNTNFP
    if cri == 'list':
        TP, FN, TN, FP = map(lambda x: len(x), [TP, FN, TN, FP])
    if not cri in ['num', 'list']: raise ValueError('Wrong cri! [num / list]')
    try: return 2*TP/float(2*TP+FP+FN)
    except: return np.nan






PATH = "../../result/drug_result/"

drug_info = pd.read_table("../../data/drug/drug_info.txt",index_col="Standard Drug Name")
drug_info = drug_info.dropna(subset=['ID'])
drug_info["Cancer(s)"] = drug_info["Cancer(s)"].replace(np.nan, "-")

cancer_list = list(pd.read_csv("../../data/cancer_list.csv")["Type"])
fda_drug = pd.read_csv("../../data/drug/fda_cancer_drug.csv", index_col="Type")
dr_set = drug_info.loc[drug_info["Status"]=="APPROVED"]
dr_set_target = filter(lambda x: "Oncology" not in drug_info.loc[x, "Pathology"].split(" | "), dr_set.index)
dr_set = dr_set.loc[dr_set_target]
cut = 0.05
result = defaultdict(list)
F1_ = {}
for ca in cancer_list:
    F1 = []
    #prms = pd.read_table("/home/kwang/PROJECT/Cancer_Driver/result/Figure5/rms/significance/%s" % ca, index_col="Drug")
    #rms = pd.read_table("/home/kwang/PROJECT/Cancer_Driver/result/Figure5/rms/%s" % ca, index_col="Drug")
    network = os.listdir("../../result/drug_result/")
    #nontar=[]
    #if ca in fda_drug.index:
    #    nontar = fda_drug.loc[ca, "drug"].split("|")
    for net in network:
        prms = pd.read_csv("../../result/drug_result/" +net+"/drug_tc_score_fdr.csv", index_col="Drug")
        #nontar = rms.loc[rms[net]==-1].index
        predict_X = []
        for val in prms[ca]:
            if val<=cut:
                predict_X.append(1)
            else: predict_X.append(0)
        predict_X = pd.Series(predict_X, index=prms.index)
        predict_X = predict_X.loc[dr_set.index]
        #predict_X = predict_X.loc[~predict_X.index.isin(nontar)]
        Y = dr_set.loc[predict_X.index, "Cancer(s)"].replace("-", 0)
        Y = Y.replace("CLINICAL CANCER", 1)
        TPFNTNFP = get_TPFNTNFP_from_strip(Y, predict_X)
        F1.append(F1_score(TPFNTNFP))
    F1_[ca] = F1
re = pd.DataFrame(F1_, index=network)
re = re.T
re.index.name = "Type"
re.to_csv("../../result/repurposing/F1_performance.csv")

