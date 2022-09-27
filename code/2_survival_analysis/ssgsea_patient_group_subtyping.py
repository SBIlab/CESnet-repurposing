import pandas as pd
import os

PATH = "../../data/Reactome/"
TYPE = os.listdir(PATH)
#TYPE.remove("group")

try:os.mkdir(PATH + "group/")
except:pass
for ty in TYPE:
    print ty
    df = pd.read_csv(PATH+ty, index_col="Term|NES")
    df = df.T
    A = []
    if len(df.index)<100:
        continue
    for pa in df.columns:
        nes = df[pa]
        g0 = nes[nes < nes.quantile(0.5)].index
        g1 = nes[nes >= nes.quantile(0.5)].index
        A.append([pa,  "|".join(g0), "|".join(g1)])
    A = pd.DataFrame(A, columns=["Term", "group1", "group2"])
    A.to_csv(PATH + "group/"+ty)
