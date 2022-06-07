# Hydrogen_Phases
A python package to compute the fractions of the various phases of Hydrogen in hydrodynamical simulations.

## Requirements

MassTensor requires:

* Python version 3.7 or later
* Numpy version 1.21.5 or later

## Optional Requirements

* h5py version 3.6.0 or later
* Matplotlib version 3.5.1 or later 

## Install

To use this package, please clone this github repository as:

`git clone https://github.com/Alex-Hill94/HydrogenPhases.git`

and add the folder to your python path. 

## Use

This package computes the relative fraction of the various phases of Hydrogen present within given gas particles of a hydrodynamical simulation. It follows the methodology presented in Blitz and Rosolowsky (2006), and Rahmati et al. (2013) for the computation of the H2 and neutral fractions respectively. From these, the fractions of (HI, HII, H2) are known. 

The computation is undertaken by the class Compute_H_Fractions within **Compute_Fractions.py**. An example can be found within **computation.py**, as well as some plotting code to illustrate the performance of the algorithm.

## Acknowledgements 

If this code is used in any published work, please acknowledge and cite the original papers, as well as this repository:

>> @ARTICLE{2006ApJ...650..933B,
       author = {{Blitz}, Leo and {Rosolowsky}, Erik},
        title = "{The Role of Pressure in GMC Formation II: The H$_{2}$-Pressure Relation}",
      journal = {\apj},
     keywords = {Galaxies: ISM, ISM: Clouds, ISM: Evolution, ISM: Molecules, Astrophysics},
         year = 2006,
        month = oct,
       volume = {650},
       number = {2},
        pages = {933-944},
          doi = {10.1086/505417},
archivePrefix = {arXiv},
       eprint = {astro-ph/0605035},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2006ApJ...650..933B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

>> @ARTICLE{2013MNRAS.430.2427R,
       author = {{Rahmati}, Alireza and {Pawlik}, Andreas H. and {Rai{\v{c}}evi{\'c}}, Milan and {Schaye}, Joop},
        title = "{On the evolution of the H I column density distribution in cosmological simulations}",
      journal = {\mnras},
     keywords = {radiative transfer, methods: numerical, galaxies: evolution, galaxies: formation, galaxies: high-redshift, intergalactic medium, Astrophysics - Cosmology and Extragalactic Astrophysics},
         year = 2013,
        month = apr,
       volume = {430},
       number = {3},
        pages = {2427-2445},
          doi = {10.1093/mnras/stt066},
archivePrefix = {arXiv},
       eprint = {1210.7808},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2427R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}





