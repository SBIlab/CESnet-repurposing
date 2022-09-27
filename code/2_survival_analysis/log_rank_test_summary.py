import pandas as pd
import numpy as np
import os

PATH = "../../result/patient_stratification/logrank/"
TYPE = os.listdir(PATH)
A = []
for ty in TYPE:
    df = pd.read_csv(PATH+ty, index_col="Term")
    A.append(df["abs_log_fold"])

re = pd.concat(A, axis=1, sort=True)
re.columns= map(lambda x: x[:-4], TYPE)
re.index.name= "Term"
re.to_csv(PATH+"ssgsea_foldchange_summary.csv")

