# CESnet-repurposing
## Description
Source codes for generating results of "A Co-essentiality Network of Cancer Driver Genes Better Prioritizes Anticancer Drugs".


## Requirements
### python2
- python (2.7.13)
- pandas (0.23.4)
- numpy (1.15.4)
- scipy (1.1.0)
- scikit-learn (0.20.4)
- networkx (2.1)
- gprofiler-official (1.0.0)
- matplotlib (2.2.5)
- seaborn (0.9.0)
- lifelines (0.19.5)

### python3
- python (3.7.9)
- pandas (1.1.3)
- numpy (1.19.2)
- scipy (1.5.2)
- gseapy (0.10.2)

## Installation
All python packages can be installed via pip (https://pypi.org/project/pip/).
pip install [package name]
e.g. pip install pandas

Gnerally, a couple of minutes is needed for installing each package.

## CESnet-repurposing
### Codes for results reproduction
Most of the codes was written in python2, but the codes related to gene set enrichment analysis(gsea) was written in python3.
Code for reproducing drug target prioritization and drug rerpurposing is provided under the './code/3_drug_candidate_analysis' folder.
Download network file, network.zip from figshare repository linked to submission page, and unzip this file under the './data/network' folder 
