import pandas as pd

df = pd.read_csv("../../data/cancer_drug_all_target_new.csv", index_col="Type")
PATH = "../../result/drug_result"
drug_info = pd.read_hdf("../../data/fitness_pandrugs_all_parsing_new.h5")
repo = drug_info.loc[drug_info["Status"] =="APPROVED",].index
d_mat = pd.read_csv(PATH + "coessentiality/drug_cancer_tc_score.csv",index_col="Drug")


cut = 10
cancer_node = {}
disease_node = {}

di_list = map(lambda x: x.split(" | "),list(drug_info.loc[repo, "Pathology"]))
disease = set(reduce(lambda x,y: x+y,di_list))

check = []
TYPE = ["DLBC", "LAML", "PAAD", "OV"]
ctype = list(d_mat.columns)
for ty in TYPE:
    ctype.remove(ty)
nodes = pd.Series(ctype + list(disease))
nodes = pd.DataFrame(nodes,columns=["names"])
nodes["index"] = nodes.index
links = []
print "nodes"
print nodes
drug_num_check = []


re_table = []
for i, ca in enumerate(ctype):
    min_table = pd.read_table("../../result/repurposing/"+ca+".txt", index_col="Drug")
    min_table = min_table.loc[repo,]
    sub_re = {}
    sub_link = []
    target = min_table.loc[min_table["coessentiality"]==-1,].index
    cand = min_table.loc[(min_table["coessentiality"] <= cut) & (min_table["coessentiality"] > 0),].index
    for cc in cand:
        re_table.append([cc, ca, drug_info.loc[cc, "Pathology"], min_table.loc[cc, "min_check"]])
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
re_table.to_csv("../../result/drug_result/repurposing_table/repurposing_list.txt", sep="\t", index=False)
#print links
links.to_csv("../../result/repurposing/repurposing_link.csv")
nodes.to_csv("../../result/repurposing/repurposing_node.csv")
