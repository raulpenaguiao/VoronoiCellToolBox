# VoronoiCellToolBox/__init__.py
r"""
VoronoiCellToolBox: Voronoi cells, secondary cones and exact second moments of lattices.

Public API (snake_case):

- :func:`voronoi_cell`, :func:`relevant_vectors`, :func:`is_generic`, :func:`reduce_metric`
- :func:`secondary_cone`, :func:`rayify`, :func:`matrixify`
- :func:`pulling_triangulation`, :func:`second_moment`, :func:`normalized_second_moment`
- :func:`chamber_second_moment`, :func:`chamber_ring`, :func:`metric_coordinates`
- :func:`plot_voronoi_cell`

The older names (``VCell``, ``normalizedChamberSecondMomentPolynomial``, ``Delsets``, ...) are
kept for backward compatibility.
"""
from VoronoiCellToolBox.voronoi_cell import (
    voronoi_cell,
    relevant_vectors,
    is_generic,
    reduce_metric,
    secondary_cone,
    rayify,
    matrixify,
    pulling_triangulation,
    second_moment,
    normalized_second_moment,
    VCell,
    eq_of_wall,
    Qform,
    Delsets,
    smMatrixSimplex,
)
from VoronoiCellToolBox.orchestration import (
    chamber_second_moment,
    chamber_ring,
    metric_coordinates,
    normalizedChamberSecondMomentPolynomial,
)
from VoronoiCellToolBox.macaulay_parsing import FormatPullingTrigMatrix
from VoronoiCellToolBox.plotting import plot_voronoi_cell

__version__ = "2.0.0"

__all__ = [
    "voronoi_cell",
    "relevant_vectors",
    "is_generic",
    "reduce_metric",
    "secondary_cone",
    "rayify",
    "matrixify",
    "pulling_triangulation",
    "second_moment",
    "normalized_second_moment",
    "chamber_second_moment",
    "chamber_ring",
    "metric_coordinates",
    "plot_voronoi_cell",
    # legacy names
    "VCell",
    "eq_of_wall",
    "Qform",
    "Delsets",
    "FormatPullingTrigMatrix",
    "smMatrixSimplex",
    "normalizedChamberSecondMomentPolynomial",
]
