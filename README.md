# Simple quasistationary method for simulations of epidemic processes with localized states: Reactivation per activity time (RAT)

This code is part of the article "[]".

[![license](https://img.shields.io/badge/licence-GPLv3-brightgreen.svg)](http://choosealicense.com/licenses/gpl-3.0/)
[![language](https://img.shields.io/badge/built%20with-Fortran-blue.svg)](https://gcc.gnu.org/fortran/)

### Fortran implementation

## Versions

* RAT - Fortran 90 (.f90)
* RAT - Python 3.7 (.py)

## Citation

Full bibliographic details: 

DOI information: 

```
@article{
}
```

## Synopsis

This code is a implementation of the RAT (Reactivation per Activity Time) quasistationary method, detailed in our [paper]. It receives a network file, containing a list of edges, as input. Some parameters can be set in the file *dyn_par.dat* (Fortran) or *dyn_par.py* (Python). 

## Inputs

You need to provide a file containing the list of edges. ID of the vertices must be enumerated sequentially as `1, 2, 3,..., N`, where `N` is the total number of vertices of the network. Here, we assume  __undirected__ and __unweighted__ networks without multiple neither self connections.

Examples of datasets at /networks. Below, you'll find a brief description of each (all networks have 10<sup>4</sup> nodes):

* *gm23.dat*: Power law degree distribution with exponent 2.3
* *gm27.dat*: Power law degree distribution with exponent 2.7
* *rrnhub.dat*: Random regular network with **m = 4** and a hub with degree = **100**

On file *dyn_par.dat*, you'll find 5 columns. Each column regards one dynamical parameter:
1. Infection rate
2. Averaging time
3. Relaxation time
4. Network file name
5. Initial seed for pseudo-random number generator only Fortran

If going for Python version, the file *dyn_par.py* will have four rows:
1. Network file
2. Infection rate
3. Relaxation time
4. Averaging time

## Outputs

For a given set of parameters, the program outputs the quasistationary distribution P(n) in **pn.dat** and three quasistationary quantities on terminal or console: density of infected vertices, dynamical susceptibility and lifespam. Details about these measures may be found in the [paper].

## Compiling and Executing
* Fortran (Gfortran): 
  * gfortran rta.f90 -o name_exec (faster, no debugging, may give multiple warnings regarding identation (insert flag -w to get rid of them))
  * gfortran rta.f90 -fcheck=all -fcheck=bounds -o name_exec (slower, for debugging)
* Fortran (Intel): 
  * ifort rta.f90 -o name_exec (faster, no debugging)
  * ifort rta.f90 -traceback -check all -o name_exec (slower, for debugging)
* Python: python rta.py

* Executing (Fortran only): ./name_exec < dyn_par.dat

Obs: The file containing the network input must be in the same folder as the .f90 or .py code.

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).
