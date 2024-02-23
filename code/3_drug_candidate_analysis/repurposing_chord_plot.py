##This code is available in python3##

import pandas as pd
import holoviews as hv
from holoviews import opts, dim
import matplotlib.pyplot as plt

hv.extension('bokeh')
link = pd.read_csv("../../result/repurposing/repurposing_link.csv")
node = pd.read_csv("../../result/repurposing/repurposing_node.csv")
del link["Unnamed: 0"]
del node["Unnamed: 0"]

tnode = node.loc[node["index"].isin(list(link["source"].unique())+list(link["target"].unique()))]

nodes = hv.Dataset(tnode, 'index')
links = link
hv.extension('matplotlib')
hv.output(fig='svg', size=250)
chord = hv.Chord((links, nodes))


chord.opts(opts.Chord(cmap='Category20', edge_cmap='Category20', edge_color=dim('source').str(), labels='names', node_color=dim('index').str()))

plt.rcParams['svg.fonttype'] = 'none'
hv.save(chord, '../../result/repurposing/repurposing_chord.svg', fmt='svg')
hv.save(chord, '../../result/repurposing/repurposing_chord.pdf', fmt='pdf')



