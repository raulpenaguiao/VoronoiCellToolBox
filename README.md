# VoronoiCellToolBox
Computes Voronoi cells as well as several useful properties and characteristis like the second moment



## Installation

This package requires a functional installation of
[SageMath](https://sagemath.org) and [Macaulay2](https://macaulay2.com). Assuming that `sage` runs Sagemath
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
`VoronoiCellToolBox` directory to be found. The `PYTHONPATH` environment variable
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

**WSL users:** the terminal link will not work directly. Instead, find your WSL IP address by running `hostname -I` and use that IP in place of `localhost`. For example, if your IP is `172.31.13.253`, open `http://172.31.13.253:8888/?token=<your-token>` in your Windows browser (copy the token from the terminal output). You will also need to start Jupyter with `--ip=0.0.0.0` to allow connections from outside WSL:
```
conda activate sage
jupyter notebook --ip=0.0.0.0 --no-browser
```

```sage
sage: from VoronoiCellToolBox.voronoi_cell import VCell
sage: VCell([[2, -1], [-1, 2]], range=2)
```

## Unit tests

```
sage -python -m unittest discover -s tests -v
```
