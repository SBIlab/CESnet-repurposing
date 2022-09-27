import networkx as nx
import pandas as pd
import numpy as np
import os, sys
import random
import pickle



alpha, gamma, num = sys.argv[1:]
#fit.txt, 0.5, 1.0, BRCA, 4
alpha = float(alpha)
gamma = float(gamma)
num = int(num)

#################
#prior knowledge
#################


with open("../../data/ukin/cgc.txt", "rb") as f:
    cgc = map(lambda x: x.strip(),f.readlines())
cgc=set(cgc)


try:
    os.mkdir("../../data/ukin/random_set/cgc_cgc/a%sr%s/" %(alpha, gamma))
except:
    pass

for nn in range(num):
    hidden = random.sample(cgc, 400)
    left = cgc - set(hidden)
    prior = list(random.sample(left, 20))
    answer = list(hidden)
    random_input = {"answer": answer, "prior": prior, "left": left}
    with open("../../data/ukin/random_set/cgc_cgc/a%sr%s/random_set_%s.pickle"% (alpha, gamma, nn), "wb") as fw:
        pickle.dump(random_input,fw,protocol=pickle.HIGHEST_PROTOCOL)








