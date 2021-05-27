# profile_binr ![](https://enktyy605xyf6r8.m.pipedream.net)

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

Once again this is a minimal example :
```python
from profile_binr import ProfileBin
import pandas as pd

# your data is assumed to contain observations as
# rows and genes as columns
data = pd.read_csv("path/to/your/data.csv")
data.head()
```
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Clec1b</th>
      <th>Kdm3a</th>
      <th>Coro2b</th>
      <th>8430408G22Rik</th>
      <th>Clec9a</th>
      <th>Phf6</th>
      <th>Usp14</th>
      <th>Tmem167b</th>
    </tr>
    <tr>
      <th>cell_id</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>HSPC_025</th>
      <td>0.0</td>
      <td>4.891604</td>
      <td>1.426148</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.599758</td>
      <td>2.954035</td>
      <td>6.357369</td>
    </tr>
    <tr>
      <th>HSPC_031</th>
      <td>0.0</td>
      <td>6.877725</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.423483</td>
      <td>1.804914</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>HSPC_037</th>
      <td>0.0</td>
      <td>0.000000</td>
      <td>6.913384</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.051659</td>
      <td>8.265465</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>LT-HSC_001</th>
      <td>0.0</td>
      <td>0.000000</td>
      <td>8.178374</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>6.419817</td>
      <td>3.453502</td>
      <td>2.579528</td>
    </tr>
    <tr>
      <th>HSPC_001</th>
      <td>0.0</td>
      <td>0.000000</td>
      <td>9.475577</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>7.733370</td>
      <td>1.478900</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>

```python
# create the binarisation instance using the dataframe
# with the index containing the cell identifier
# and the columns being the gene names
probin = ProfileBin(data)

# compute the criteria used to binarise/normalise the data :
# This method uses a parallel implementation, you can specify the 
# number of workers with an integer
probin.fit(8) # train using 8 threads

# Look at the computed criteria
probin.criteria.head(8)
```
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Dip</th>
      <th>BI</th>
      <th>Kurtosis</th>
      <th>DropOutRate</th>
      <th>MeanNZ</th>
      <th>DenPeak</th>
      <th>Amplitude</th>
      <th>Category</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Clec1b</th>
      <td>0.358107</td>
      <td>1.635698</td>
      <td>54.017736</td>
      <td>0.876208</td>
      <td>1.520978</td>
      <td>-0.007249</td>
      <td>8.852181</td>
      <td>ZeroInf</td>
    </tr>
    <tr>
      <th>Kdm3a</th>
      <td>0.000000</td>
      <td>2.407548</td>
      <td>-0.784019</td>
      <td>0.326087</td>
      <td>3.847940</td>
      <td>0.209239</td>
      <td>10.126676</td>
      <td>Bimodal</td>
    </tr>
    <tr>
      <th>Coro2b</th>
      <td>0.000000</td>
      <td>2.320060</td>
      <td>7.061604</td>
      <td>0.658213</td>
      <td>2.383819</td>
      <td>0.004597</td>
      <td>9.475577</td>
      <td>ZeroInf</td>
    </tr>
    <tr>
      <th>8430408G22Rik</th>
      <td>0.684454</td>
      <td>3.121069</td>
      <td>21.729044</td>
      <td>0.884058</td>
      <td>2.983472</td>
      <td>0.005663</td>
      <td>9.067857</td>
      <td>ZeroInf</td>
    </tr>
    <tr>
      <th>Clec9a</th>
      <td>1.000000</td>
      <td>2.081717</td>
      <td>140.089285</td>
      <td>0.965580</td>
      <td>2.280293</td>
      <td>-0.009361</td>
      <td>9.614233</td>
      <td>Discarded</td>
    </tr>
    <tr>
      <th>Phf6</th>
      <td>0.000000</td>
      <td>1.988667</td>
      <td>-1.389024</td>
      <td>0.035628</td>
      <td>5.025501</td>
      <td>2.017547</td>
      <td>10.135226</td>
      <td>Bimodal</td>
    </tr>
    <tr>
      <th>Usp14</th>
      <td>0.000000</td>
      <td>2.208080</td>
      <td>-1.224987</td>
      <td>0.007850</td>
      <td>6.109964</td>
      <td>8.245570</td>
      <td>11.088750</td>
      <td>Bimodal</td>
    </tr>
    <tr>
      <th>Tmem167b</th>
      <td>0.000000</td>
      <td>2.430813</td>
      <td>0.093023</td>
      <td>0.393720</td>
      <td>3.448331</td>
      <td>0.072982</td>
      <td>9.486826</td>
      <td>Bimodal</td>
    </tr>
  </tbody>
</table>


```python
# get binarised data (alternatively .binarise()):
my_bin = probin.binarize()
my_bin.head()
```
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Clec1b</th>
      <th>Kdm3a</th>
      <th>Coro2b</th>
      <th>8430408G22Rik</th>
      <th>Clec9a</th>
      <th>Phf6</th>
      <th>Usp14</th>
      <th>Tmem167b</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>HSPC_025</th>
      <td>NaN</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>HSPC_031</th>
      <td>NaN</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HSPC_037</th>
      <td>NaN</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>LT-HSC_001</th>
      <td>NaN</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>HSPC_001</th>
      <td>NaN</td>
      <td>0.0</td>
      <td>1.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>


```python
# idem for normalised data :
my_norm = probin.normalize()
my_norm.head()
```
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Clec1b</th>
      <th>Kdm3a</th>
      <th>Coro2b</th>
      <th>8430408G22Rik</th>
      <th>Clec9a</th>
      <th>Phf6</th>
      <th>Usp14</th>
      <th>Tmem167b</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>HSPC_025</th>
      <td>0.0</td>
      <td>9.786196e-01</td>
      <td>0.184102</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.000801</td>
      <td>8.318176e-05</td>
      <td>9.999970e-01</td>
    </tr>
    <tr>
      <th>HSPC_031</th>
      <td>0.0</td>
      <td>9.999981e-01</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.000462</td>
      <td>8.084114e-07</td>
      <td>6.874397e-11</td>
    </tr>
    <tr>
      <th>HSPC_037</th>
      <td>0.0</td>
      <td>4.408417e-09</td>
      <td>0.892449</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.000145</td>
      <td>9.999940e-01</td>
      <td>6.874397e-11</td>
    </tr>
    <tr>
      <th>LT-HSC_001</th>
      <td>0.0</td>
      <td>4.408417e-09</td>
      <td>1.000000</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.991865</td>
      <td>6.230178e-04</td>
      <td>1.599753e-04</td>
    </tr>
    <tr>
      <th>HSPC_001</th>
      <td>0.0</td>
      <td>4.408417e-09</td>
      <td>1.000000</td>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.999865</td>
      <td>2.171153e-07</td>
      <td>6.874397e-11</td>
    </tr>
  </tbody>
</table>


## References

* Béal J, Montagud A, Traynard P, Barillot E and Calzone L (2019) *Personalization of Logical Models With Multi-Omics Data Allows Clinical Stratification of Patients*. Front. Physiol. 9:1965. doi:[10.3389/fphys.2018.01965](https://doi.org/10.3389/fphys.2018.01965)
