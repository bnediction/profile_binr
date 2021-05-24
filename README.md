# profile_binr

The PROFILE methodology for the binarisation and normalisation of RNA-seq data.

This is a Python interface to a set of normalisation and binarisation functions for
RNA-seq data originally written in R.

This software package is based on the methodology developed by Beal, Jonas; Montagud, Arnau;
Traynard, Pauline; Barillot, Emmanuel; and Calzone, Laurence at [Computational Systems Biology of Cancer team at Institut Curie](https://sysbio.curie.fr/) 
([contact-sysbio@curie.fr](mailto:contact-sysbio@curie.fr)).

[This is the repository](https://github.com/sysbio-curie/PROFILE) containing the original implementation in
Rmarkdown notebooks.

# Installation

This software has only been tested in Debian-based GNU/Linux distributions, it should in principle
work on any *nix system.

## Prerequisites (non-python deoendencies)

* R, version 4.0.2 (2020-06-22) -- "Taking Off Again"
* git (temporarily, to be removed)
* make
* g++
* gfortran


## Using pip

This is a barebones functional example. We recommend installing within a Python virtual environment.

```bash
pip install git+https://github.com/bnediction/profile_binr
```

# Usage

Once again this is a minimal example :
```python
from profile_binr.wrappers.probinr import ProfileBin
import pandas as pd

# your data is assumed to contain observations as
# rows and genes as columns
data = pd.read_csv("path/to/your/data.csv")

probin = ProfileBin(data)

# compute the criteria used to binarise/normalise the data :
probin.fit()

# get binarised data :
my_bin = probin.binarize()
# alternatively :
# my_bin = probin.binarise()

# idem for normalised data :
my_norm = probin.normalize()
```

# References 

Please use the following bibtex entries to cite the original author's work :


```latex
@article{Beal2019,
abstract = {Logical models of cancer pathways are typically built by mining the literature for relevant experimental observations. They are usually generic as they apply for large cohorts of individuals. As a consequence, they generally do not capture the heterogeneity of patient tumors and their therapeutic responses. We present here a novel framework, referred to as PROFILE, to tailor logical models to a particular biological sample such as a patient tumor. This methodology permits to compare the model simulations to individual clinical data, i.e., survival time. Our approach focuses on integrating mutation data, copy number alterations (CNA), and expression data (transcriptomics or proteomics) to logical models. These data need first to be either binarized or set between 0 and 1, and can then be incorporated in the logical model by modifying the activity of the node, the initial conditions or the state transition rates. The use of MaBoSS, a tool based on Monte-Carlo kinetic algorithm to perform stochastic simulations on logical models results in model state probabilities, and allows for a semi-quantitative study of the model phenotypes and perturbations. As a proof of concept, we use a published generic model of cancer signaling pathways and molecular data from METABRIC breast cancer patients. For this example, we test several combinations of data incorporation and discuss that, with these data, the most comprehensive patient-specific cancer models are obtained by modifying the nodes' activity of the model with mutations, in combination or not with CNA data, and altering the transition rates with RNA expression. We conclude that these model simulations show good correlation with clinical data such as patients' Nottingham prognostic index (NPI) subgrouping and survival time. We observe that two highly relevant cancer phenotypes derived from personalized models, Proliferation and Apoptosis, are biologically consistent prognostic factors: patients with both high proliferation and low apoptosis have the worst survival rate, and conversely. Our approach aims to combine the mechanistic insights of logical modeling with multi-omics data integration to provide patient-relevant models. This work leads to the use of logical modeling for precision medicine and will eventually facilitate the choice of patient-specific drug treatments by physicians.},
author = {Beal, Jonas and Montagud, Arnau and Traynard, Pauline and Barillot, Emmanuel and Calzone, Laurence},
doi = {10.3389/fphys.2018.01965},
file = {:home/gml/Documents/Master/M1_BIBS/S2/TER/biblio/personal-models-omics.pdf:pdf},
issn = {1664042X},
journal = {Frontiers in Physiology},
keywords = {Breast cancer,Data discretization,Logical models,Personalized mechanistic models,Personalized medicine,Stochastic simulations},
number = {JAN},
title = {{Personalization of logical models with multi-omics data allows clinical stratification of patients}},
volume = {10},
year = {2019}
}
@article{Beal2019a,
abstract = {Logical models of cancer pathways are typically built by mining the literature for relevant experimental observations. They are usually generic as they apply for large cohorts of individuals. As a consequence, they generally do not capture the heterogeneity of patient tumors and their therapeutic responses. We present here a novel framework, referred to as PROFILE, to tailor logical models to a particular biological sample such as a patient tumor. This methodology permits to compare the model simulations to individual clinical data, i.e., survival time. Our approach focuses on integrating mutation data, copy number alterations (CNA), and expression data (transcriptomics or proteomics) to logical models. These data need first to be either binarized or set between 0 and 1, and can then be incorporated in the logical model by modifying the activity of the node, the initial conditions or the state transition rates. The use of MaBoSS, a tool based on Monte-Carlo kinetic algorithm to perform stochastic simulations on logical models results in model state probabilities, and allows for a semi-quantitative study of the model phenotypes and perturbations. As a proof of concept, we use a published generic model of cancer signaling pathways and molecular data from METABRIC breast cancer patients. For this example, we test several combinations of data incorporation and discuss that, with these data, the most comprehensive patient-specific cancer models are obtained by modifying the nodes' activity of the model with mutations, in combination or not with CNA data, and altering the transition rates with RNA expression. We conclude that these model simulations show good correlation with clinical data such as patients' Nottingham prognostic index (NPI) subgrouping and survival time. We observe that two highly relevant cancer phenotypes derived from personalized models, Proliferation and Apoptosis, are biologically consistent prognostic factors: patients with both high proliferation and low apoptosis have the worst survival rate, and conversely. Our approach aims to combine the mechanistic insights of logical modeling with multi-omics data integration to provide patient-relevant models. This work leads to the use of logical modeling for precision medicine and will eventually facilitate the choice of patient-specific drug treatments by physicians.},
author = {Beal, Jonas and Montagud, Arnau and Traynard, Pauline and Barillot, Emmanuel and Calzone, Laurence},
doi = {10.3389/fphys.2018.01965},
file = {:home/gml/Documents/Master/M1_BIBS/S2/TER/biblio/additional-personal-models-omics.pdf:pdf},
issn = {1664042X},
journal = {Frontiers in Physiology},
keywords = {Breast cancer,Data discretization,Logical models,Personalized mechanistic models,Personalized medicine,Stochastic simulations},
number = {JAN},
pages = {1--23},
title = {{Personalization of logical models with multi-omics data allows clinical stratification of patients}},
volume = {10},
year = {2019}
}

```


