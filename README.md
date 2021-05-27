# profile_binr 

![](https://enktyy605xyf6r8.m.pipedream.net)
[![Join the chat at https://gitter.im/bnediction/profile_binr](https://badges.gitter.im/bnediction/profile_binr.svg)](https://gitter.im/bnediction/profile_binr?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

The PROFILE methodology for the binarisation and normalisation of RNA-seq data.

This is a Python interface to a set of normalisation and binarisation functions for
RNA-seq data originally written in R.

This software package is based on the methodology developed by Beal, Jonas; Montagud, Arnau;
Traynard, Pauline; Barillot, Emmanuel; and Calzone, Laurence at [Computational Systems Biology of Cancer team at Institut Curie](https://sysbio.curie.fr/)
([contact-sysbio@curie.fr](mailto:contact-sysbio@curie.fr)).
It generalizes and offers a Python interface of the original implementation in Rmarkdown notebooks available at https://github.com/sysbio-curie/PROFILE.

## Installation

### Using conda

The tool can be installed using the Conda package [profile_binr](https://anaconda.org/colomoto/profile_binr) in the `colomoto` channel. Note that some of its dependencies requires the `conda-forge` channel.

```
conda install -c conda-forge colomoto::profile_binr
```

### Using pip

#### Requirements

* R (≥4.0)
* R packages:
    * mclust
    * diptest
    * moments
    * magrittr
    * tidyr
    * dplyr
    * tibble
    * bigmemory
    * doSNOW
    * foreach
    * glue


```
pip install profile_binr
```

## Usage

A minimal example of the binarization suite. Take a look at `notebooks/` for more details.
```python
from profile_binr import ProfileBin
import pandas as pd

# your data is assumed to contain observations as
# rows and genes as columns
data = pd.read_csv("path/to/your/data.csv")
data.head()

# create the binarisation instance using the dataframe
# with the index containing the cell identifier
# and the columns being the gene names
probin = ProfileBin(data)

# compute the criteria used to binarise/normalise the data :
# This method uses a parallel implementation, you can specify the 
# number of workers with an integer
probin.fit(8) # train using 8 threads

# Look at the computed criteria
probin.criteria

# get binarised data (alternatively .binarise()):
my_bin = probin.binarize()
my_bin.head()

# idem for normalised data :
my_norm = probin.normalize()
my_norm.head()
```

## References

* Béal J, Montagud A, Traynard P, Barillot E and Calzone L (2019) *Personalization of Logical Models With Multi-Omics Data Allows Clinical Stratification of Patients*. Front. Physiol. 9:1965. doi:[10.3389/fphys.2018.01965](https://doi.org/10.3389/fphys.2018.01965)
