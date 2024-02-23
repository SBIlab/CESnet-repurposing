import pandas as pd
import numpy as np
import os
import lifelines as lfl



PATH = PATH = "../../data/Reactome/group/"
TYPE = os.listdir(PATH)
for ty in TYPE:
    print ty
    gr = pd.read_csv(PATH + ty, index_col="Term")
    del gr["Unnamed: 0"]


    pat = pd.read_csv("../../data/TCGA_patient/TCGA-" +ty[:-4]+ "/clinical_patient_table.csv",
                      skiprows=[0,2], index_col="bcr_patient_barcode")[["vital_status", "days_to_last_followup", "days_to_death"]]
    pat = pat.replace(["[Not Available]","[Not Applicable]",],"0")
    pat = pat.replace(["[Completed]",],"-1")
    pat = pat.replace("[Discrepancy]", "-9999999") # drop negative value
    pat["T"] = pat["days_to_death"].astype(int) + pat["days_to_last_followup"].astype(int)
    pat["E"] = 0
    pat.loc[pat["vital_status"]=="Dead","E"] = 1
    pat = pat.loc[pat["T"]>=0,:]
    result = []

    for pa in gr.index:
        g0 = map(lambda x: x[:12], gr.loc[pa,"group1"].split("|"))
        g1 = map(lambda x: x[:12], gr.loc[pa,"group2"].split("|"))
        group = []
        for pp in pat.index:
            if pp in g0:
                group.append(0)
            elif pp in g1:
                group.append(1)
            else:
                group.append(2)
        pat["group"] = group
        suv = pat.loc[pat["group"] != 2, ["T", "E", "group"]]
        try:
            groups = suv["group"]
            ix = (groups == 0)
            re = lfl.statistics.logrank_test(suv.loc[~ix, "T"], suv.loc[ix, "T"], event_observed_A=suv.loc[~ix, "E"], event_observed_B=suv.loc[ix, "E"])
            result.append([pa, re.p_value, re.test_statistic])
            #ind.append(pa)
        except:
            result.append([pa, np.nan, np.nan])
        del pat["group"]
    result = pd.DataFrame(result, columns=["Term", "p-value", "statistics"])
    result = result.set_index("Term")
    result.to_csv("../../result/patient_stratification/logrank/"+ty)





