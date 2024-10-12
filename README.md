# fitness-landscapes

Repository for 'Impact of population size on early adaptation in rugged fitness landscapes' (2023) by R. Servajean and A.-F. Bitbol available here: https://royalsocietypublishing.org/doi/full/10.1098/rstb.2022.0045
and 'Impact of spatial structure on early and long-term adaptation in rugged fitness landscapes' (2024) by R. Servajean, A. Alexandre and A.-F. Bitbol available here: https://www.biorxiv.org/content/10.1101/2024.09.23.614481.abstract

## Getting started ##

Clone this repository on your local machine by running:

```bash
git clone git@github.com:Bitbol-Lab/fitness-landscapes.git
``` 
 

Executing the following line runs a working example of the Moran and Wright-Fisher walks:
```bash
python moran_and_wright_fisher_walks.py
```

Executing the following line runs a working example of the Moran and Wright-Fisher star walks:
```bash
python star_walk.py
``` 

## Requirements ##

In order to use the function `moran_wf_walk` in `moran_and_wright_fisher_walks.py`, or `star_walk` in `star_walk.py`, NumPy is required.


## Usage ##

In the file `moran_and_wright_fisher_walks.py`,
`
moran_wf_walk
`
runs either Moran or Wright-Fisher walks in a given fitness landscape.

In the file `star_walk.py`,
`
star_walk
`
runs either Moran or Wright-Fisher star walks in a given fitness landscape.
