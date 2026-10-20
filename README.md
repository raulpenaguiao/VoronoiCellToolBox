# VoronoiCellToolBox

Exact Voronoi cells, secondary cones and second moments of lattices, in [SageMath](https://sagemath.org) and [Macaulay2](https://macaulay2.com).

A positive definite symmetric matrix `Q` (a Gram matrix) defines the metric `x^T Q y` on `R^d`, and with it a lattice `Z^d`. Given `Q`, the package computes, exactly over the rationals:

- the **Voronoi cell** `V_Q`, its relevant vectors, and whether `Q` is generic (the cell is primitive). Every cell is certified to have volume 1, so no relevant vector is missed;
- the **secondary cone** of `Q`: all metrics whose Voronoi cell has the same combinatorics, as inequalities and rays;
- a **pulling triangulation** and the **second moment** `SM(Q) = ∫_{V_Q} x^T Q x dx`, with the dimensionless `G(Q) = SM(Q) / (d det(Q)^{1/d})` of the lattice quantizer problem;
- the **chamber second moment**: `SM` as a rational function of the entries of the metric, valid on the whole secondary cone. The symbolic integration runs in Macaulay2.

Metrics in any chamber are accepted: they are LLL-reduced internally and the results are mapped back.

## Installation

Requirements: SageMath (tested with 10.6 and 10.10) and Macaulay2 (tested with 1.22), with the `M2` executable on the `PATH`. Then

    sage -pip install git+https://github.com/raulpenaguiao/VoronoiCellToolBox

or, from a clone of this repository, `make install` (which runs `sage -pip install --upgrade .`).

## Quick start

```sage
sage: from VoronoiCellToolBox import *
sage: Q = [[2, -1], [-1, 2]]                      # the hexagonal lattice
sage: second_moment(Q), normalized_second_moment(Q)
(5/18, 5/108*sqrt(3))
sage: [matrixify(list(r), 2) for r in rayify(secondary_cone(Q))]
[[[0, 0], [0, 1]], [[1, -1], [-1, 1]], [[1, 0], [0, 0]]]
sage: f = chamber_second_moment(Q); f            # variables: q_0 = q11, q_1 = q12, q_2 = q22
(2*q_0*q_1^2 + 2*q_1^3 - q_0^2*q_2 + 2*q_1^2*q_2 - q_0*q_2^2)/(12*q_1^2 - 12*q_0*q_2)
sage: f(*metric_coordinates([[3, -1], [-1, 5]])) == second_moment([[3, -1], [-1, 5]])
True
```

The notebook [`demo.ipynb`](demo.ipynb) goes from dimension 2 to the three chambers in dimension 4 (open it with the SageMath kernel, e.g. `sage -n jupyter`).

## Main functions

| Function | Returns |
|---|---|
| `voronoi_cell(Q)` | the Voronoi cell, a Sage `Polyhedron` |
| `relevant_vectors(Q)`, `is_generic(Q)` | the facet vectors; whether the cell is primitive |
| `secondary_cone(Q)`, `rayify(C)`, `matrixify(v, d)` | the secondary cone (coordinates: upper triangle of `Q`, row by row); its rays; a ray as a matrix |
| `pulling_triangulation(P)` | a triangulation of a polytope without new vertices |
| `second_moment(Q)`, `normalized_second_moment(Q)` | `SM(Q)` (rational) and `G(Q)` (exact, symbolic) |
| `chamber_second_moment(Q)` | `SM` on the secondary cone of `Q`, in `chamber_ring(d)`; evaluate with `f(*metric_coordinates(Q'))` |
| `plot_voronoi_cell(Q)` | a picture in dimension 2 or 3 |
| `reduce_metric(Q)` | the LLL-reduced metric `U^T Q U` and `U` |

Every function has a docstring with examples. Older names (`VCell`, `normalizedChamberSecondMomentPolynomial`, ...) still work.

## Tests

    make test       # unit and regression tests: sage -python -m pytest tests -v (about a minute)
    make doctest    # the docstring examples: sage -t --long VoronoiCellToolBox/

The regression tests compare against an independent implementation of the second moment (`tests/independent_checker.py`), against the secondary cones of the three Delone triangulations in dimension 4 from Zimmermann's thesis, and against the known values `G(A_4^*) ≈ 0.077559` and `G(D_4) ≈ 0.0766032346`.

## Authors and license

Yelena Mandelshtam and Raul Penaguiao. Released under the GNU General Public License, version 3 or later (see `LICENSE`).
