"""
The PROFILE methodology for the binarisation and normalisation of RNA-seq data.
"""

from .wrappers.probinr import ProfileBin
from .utils.normalization import log_transform, normalize, log_normalize

__version__ = "0.1.3a"
__author__ = "Gustavo Magaña López"
__credits__ = "BNediction ; Institut Curie"
