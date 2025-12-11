# Time Domain F-statistics (TDFstat) pipeline

This repository contains Time Domain F-statistics (TDFstat) pipeline used to search for almost-monochromatic, continuous gravitational waves in data from interferometric detectors.
Supported detectors are H1 (Hanford), L1 (Livingston) of LIGO and V1 (Cascina) of Virgo.
It is authored by members of the [Polgraw group](https://polgraw.camk.edu.pl/) within the [Virgo collaboration](https://www.virgo-gw.eu/).

See the [documentation page](https://polgraw.github.io/TDFstat/) for more details.

> :warning: Repository under construction! We are moving here from [the old repository](https://github.com/mbejger/polgraw-allsky).

The repository contains a submodule (`gridgen`). To clone with submodule use:
```
git clone --recurse-submodules https://github.com/Polgraw/TDFstat.git
#or 
git clone https://github.com/Polgraw/TDFstat.git
git submodule init
git submodule update
```

## Contributors

In alphabetic order: 

* Pia Astone 
* Michał Bejger
* Jan Bolek
* Paweł Ciecieląg
* Orest Dorosh
* Aleksander Garus
* Andrzej Królak
* Máté Ferenc Nagy-Egri
* Maciej Piętka
* Andrzej Pisarski 
* Gevorg Poghosyan
* Magdalena Sieniawska 
* Rafał Skrzypiec


This repository gathers several elements of the pipeline: 

* Extraction of narrow band short time series from SFDB: `extract_band_hdf_td`
* Creation of long duration time series: `genseg`
* Parameter space grid generation:`gridgen` (this is a copy of external repository https://github.com/apisarski/gridgen)
* Coherent search for candidate signals in narrowband time segments: `search`

TO BE ADDED LATER

* Search for coincidences among different time segments: `coincidences` 
* Estimation of false alarm probability of coincidences: `fap`
* Followup of interesting outliers: `followup`
* Sensitivity upper limits: in `sensitivity-scripts`
* [Test data](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz).  

We also include some external libraries in the utils directory:

* [iniparser](https://github.com/ndevilla/iniparser) by N. Devillard & E. Leblond,
* [Sleef vector library](https://sleef.org/) by N. Shibata et al.

Many thanks to the authors!
