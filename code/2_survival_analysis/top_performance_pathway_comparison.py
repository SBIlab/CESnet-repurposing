import pandas as pd
import numpy as np
import os


PATH = "../../data/network/"
network = os.listdir(PATH)
lo = pd.read_csv("../../result/patient_stratification/logrank/ssgsea_log_rank_summary.csv", index_col="Term")


print network
for net in network:
    print net
    TYPE = os.listdir(PATH + net +"/Reactome/")
    ### check point ###
    TYPE = filter(lambda x: x[:-9] in lo.columns,TYPE)
    result = {}
    re_index = ["fdr_0.001"]
    for ty in TYPE:
        gsea = pd.read_table(PATH + net +"/Reactome/"+ty, index_col="Term")
        gsea = gsea.loc[gsea["fdr"]<=0.001,]
        gsea = gsea.loc[gsea["es"]>0.,]
        gsea = gsea.loc[gsea["nes"]>0.,]
        result_list = []
        ### check point ###
        if len(lo.loc[gsea.index, ty[:-9]].sort_values().index)>0:
            result_list.append(lo.loc[gsea.index, ty[:-9]].min())
        else: result_list.append(np.nan)
        result[ty[:-9]] = result_list
    result = pd.DataFrame(result, index=re_index)
    result.index.name = "top"
    print result
    try:
        os.mkdir("../../result/patient_stratification/network/" + net)
    except: pass
    result.to_csv("../../result/patient_stratification/network/" + net + "/network_patients_min_logrank_result.csv")