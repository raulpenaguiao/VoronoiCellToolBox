comma := ,

.PHONY: all install build test doctest check

all: install

install:
	sage -pip install --upgrade .

build:
	sage --python setup.py build

# Unit and regression tests (about a minute; the d=4 checks dominate)
test:
	sage -python -m pytest tests -v

# The examples in the docstrings (those tagged "optional - macaulay2" only if M2 is installed)
doctest:
	sage -t --long $(if $(shell command -v M2),--optional=sage$(comma)macaulay2) VoronoiCellToolBox/

check: test doctest
