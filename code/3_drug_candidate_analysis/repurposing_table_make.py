import pandas as pd
import numpy as np
import os

drug_info = pd.read_table("../../data/drug/drug_info.txt", index_col="Standard Drug Name")
drug_info = drug_info.dropna(subset=["ID"])
repo = drug_info.loc[drug_info["Status"] =="APPROVED",]#.index
repo = filter(lambda x: "Oncology" not in repo.loc[x, "Pathology"].split(" | "), repo.index)

cut = 0.05


di_list = map(lambda x: x.split(" | "),list(drug_info.loc[repo, "Pathology"]))
disease = set(reduce(lambda x,y: x+y,di_list))

check = []
cancer_list = list(pd.read_csv("../../data/cancer_list.csv")["Type"])
nodes = pd.Series(cancer_list + list(disease))
nodes = pd.DataFrame(nodes,columns=["names"])
nodes["index"] = nodes.index
links = []
print "nodes"
print nodes
drug_num_check = []

re_table = []

network = os.listdir("../../result/drug_result")
network.remove("co-essentiality")

min_table = pd.read_csv("../../result/drug_result/co-essentiality/drug_tc_score_fdr.csv", index_col = "Drug")
for i, ca in enumerate(cancer_list):
    rp_list = []
    for net in network:
        rp_list.append(pd.read_csv("../../result/drug_result/%s/drug_tc_score_fdr.csv" % net, index_col="Drug").loc[repo,ca])
    rpnet = pd.concat(rp_list, axis=1)
    rpnet.index = repo
    rpnet.column = network
    rpcheck = (rpnet < cut).sum(axis=1)
    min_table = min_table.loc[repo,]
    sub_re = {}
    sub_link = []
    cand = min_table.loc[min_table[ca] <= cut,].index
    for cc in cand:
        if rpcheck.loc[cc] >0:
            re_table.append([cc, ca, drug_info.loc[cc, "Pathology"], 0])
        else: re_table.append([cc, ca, drug_info.loc[cc, "Pathology"], 0])
    drug_num_check += list(cand)

    check += list(cand)
    repo_info = drug_info.loc[cand, "Pathology"]
    for j,rr in enumerate(repo_info):
        pathology = rr.split(" | ")
        pathology = filter(lambda x: x != "-" ,pathology)
        for pat in pathology:
            edge = ca + "-" + pat
            if edge in sub_re.keys():
                sub_re[edge] += 1
            else:
                sub_re[edge] = 1

    for ss in sub_re.keys():
        sl = ss.split("-")
        sub_link.append([nodes.loc[nodes["names"] == sl[0],"index"].iloc[0], nodes.loc[nodes["names"] == sl[1],"index"].iloc[0], sub_re[ss]])
    sub_link = pd.DataFrame(sub_link, columns=["source", "target", "value"])
    links.append(sub_link)
links = pd.concat(links)
print "drug_pair_num:", len(check)
print "drug_num:", len(set(drug_num_check))
re_table = pd.DataFrame(re_table, columns=["drug","repurpose","origin","min_check"])
re_table.to_csv("../../result/repurposing/repurposing_list_pval.txt", sep="\t", index=False)
links.to_csv("../../result/repurposing/repurposing_link_min_pval.csv")
nodes.to_csv("../../result/repurposing/repurposing_node_min_pval.csv")
