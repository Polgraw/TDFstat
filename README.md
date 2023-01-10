# Time Domain F-statistics (TDFstat) pipeline

This repository contains Time Domain F-statistics (TDFstat) pipeline used to search for almost-monochromatic, continuous gravitational waves in data from interferometric detectors.
It is authored by members of the [Polgraw group](https://polgraw.camk.edu.pl/) within the [Virgo collaboration](https://www.virgo-gw.eu/).

See the [documentation page](https://polgraw.github.io/TDFstat/) for more details.

> :warning: Repository under construction! We are moving here from [the old repository](https://github.com/mbejger/polgraw-allsky).

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

* Generation of the initial data: `genseg`
* Parameter space grid generation:`gridgen`
* Coherent search for candidate signals in narrowband time segments: `search`
* Search for coincidences among different time segments: `coincidences` 
* Estimation of false alarm probability of coincidences: `fap`
* Followup of interesting outliers: `followup`
* Sensitivity upper limits: in `sensitivity-scripts`
* [Test data](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz).  

This pipeline searches for continuous gravitational wave signals in time-domain data using the F-statistic on data from a network of detectors (currently the available detectors are H1 (Hanford), L1 (Livingston) of LIGO and V1 (Cascina) of Virgo)  

