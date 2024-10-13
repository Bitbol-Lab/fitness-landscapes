# fitness-landscapes

Repository for 'Impact of population size on early adaptation in rugged fitness landscapes' (2023) by R. Servajean and A.-F. Bitbol available here: https://royalsocietypublishing.org/doi/full/10.1098/rstb.2022.0045 (referred to as [SER23] below)
and 'Impact of spatial structure on early and long-term adaptation in rugged fitness landscapes' (2024) by R. Servajean, A. Alexandre and A.-F. Bitbol available here: https://www.biorxiv.org/content/10.1101/2024.09.23.614481.abstract (referred to as [SER24] below)

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

## Landscapes ##

The 200,000 LK landscapes (with L=3 and K=1) used in both [SER23] and [SER24] are stored in a file called `landscapes.pickle` available here: https://drive.google.com/drive/folders/1Er1wZb_FQr2e-NWZU279jvjZxgnWJVWc?usp=sharing

Each landscape is associated with an index (started from 0). The following line of Python code provides access to the *i*th landscape (as a dictionnary where a key is a sequence and the value is the corresponding fitness):
```
landscape_i = pickle.load(open("landscapes.pickle", "rb"))[i]
```

Landscapes labelled A and B in both [SER23] and [SER24] correspond to landscapes 9 and 3, respectively. Landscape C in [SER24] corresponds to landscape 142706.
