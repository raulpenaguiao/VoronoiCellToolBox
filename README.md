# VoronoiCellToolBox
Computes Voronoi cells as well as several useful properties and characteristis like the second moment



## Installation

This package requires a functional installation of
[SageMath](https://sagemath.org). Assuming that `sage` runs Sagemath
and that you have write-permission on the sagemath install, you should
be able to install this package into sage with something like

    sage --pip install git+https://github.com/raulpenaguiao/VoronoiCellToolBox

If you do not have write permission on the sagemath install itself, you
can get a copy of the repository by something like

    git clone https://github.com/raulpenaguiao/VoronoiCellToolBox.git
    
With a command like

    sage --python setup.py build

it is still possible to build the module. You will just have to make sure
that in Sagemath, the variable `sys.path` contains a path that allows the built
`riemann_theta` directory to be found. The `PYTHONPATH` environment variable
may be useful for this. One can also install the package as a "user" package
by executing

    sage --python setup.py install --user
    
but pay attention to the warning: sage is by default configured to not
look at per-user environments. 

## Usage

One can follow through the examples in the jupyter-notebook by navigating on a terminal to the VoronoiCellToolBox directory and running
```
sage -n jupyter
```
If the webpage displays an error `Access to the file was denied` then simply open in the browser [a new jupyter notebook](http://localhost:8888).

sage: from VoronoiCellToolBox.voronoi_cell import VCell
sage: VCell([[2, -1], [-1, 2]], range=2)


## Unit tests

```
sage -python -m unittest discover -s tests -v
```
