# Simple quasistationary method for simulations of epidemic processes with localized states: Reactivation per activity time (RAT)

This code is part of the article "[]".

[![license](https://img.shields.io/badge/licence-GPLv3-brightgreen.svg)](http://choosealicense.com/licenses/gpl-3.0/)
[![language](https://img.shields.io/badge/built%20with-Fortran-blue.svg)](https://gcc.gnu.org/fortran/)

### Fortran implementation

## Versions

[Fortran implementation]

[Python implementation]

## Citation

Full bibliographic details: 

DOI information: 

```
@article{
}
```

## Synopsis

This code is a implementation of the RAT (Reactivation per Activity Time) quasistationary method, detailed in our [paper]. It receives a network file, containing a list of edges, as input. Some dynamical parameters can be set in the main code's header. They are:

$\lambda $

## Dataset input

You need to provide a file containing the list of edges (__in__ and __out__, two collumns). ID of the vertices must be enumerated sequentially as `1, 2, 3,..., N`, where `N` is the total number of vertices of the network. Here, we assume  __undirected__ and __unweighted__ networks without multiple neither self connections.

Consider, for example, a network with `N=5` vertices represented by:

```
1,2
1,3
2,4
2,5
3,4
```

Examples of datasets and their specifications are available at http://goo.gl/Bm3VsR.

## Installation

In Linux and OSX, it is simple: just type ``make`` in the terminal in the *.f90 directory. If you need debugging, use ``make c=1``.

For Windows, however, you must compile all mod*.f90 files and the program code dynamics.f90. An example is:

```gfortran mod_read_tools.f90 mod_random.f90 mod_netdata.f90 dynamics.f90 -o dynamics```


## Use

If you want to manually input the dynamical parameters, just type:

```./dynamics <edges_file> <output_file>```

where ``<output_file>`` will be written with the average density of infected vertices versus time.

Alternatively, use (Linux):

```bash run.sh <edges_file> <output_file> <number of samples> <infection rate lambda> <maximum time steps> <fraction of infected vertices (initial condition)>```

_Example:_

```bash run.sh edges/s01.edges.dat "s01.lb0.002_100-samples.dat" 100 0.002 1000000 0.5```

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).
