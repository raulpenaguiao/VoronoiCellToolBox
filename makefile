.PHONY: all install build test doctest check

all: install

install:
	sage -pip install --upgrade .

build:
	sage --python setup.py build

# Unit and regression tests (about a minute; the d=4 checks dominate)
test:
	sage -python -m pytest tests -v

# The examples in the docstrings
doctest:
	sage -t --long VoronoiCellToolBox/

check: test doctest
