# VoronoiCellToolBox/__init__.py
r"""
VoronoiCellToolBox: Voronoi cells, secondary cones and exact second moments of lattices.

Public API (snake_case):

- :func:`voronoi_cell`, :func:`relevant_vectors`, :func:`is_generic`, :func:`reduce_metric`
- :func:`secondary_cone`, :func:`rayify`, :func:`matrixify`
- :func:`pulling_triangulation`, :func:`second_moment`, :func:`normalized_second_moment`
- :func:`chamber_second_moment`, :func:`chamber_ring`, :func:`metric_coordinates`
- :func:`plot_voronoi_cell`
- graph lattices: :func:`riemann_matrix`, :func:`cycle_basis_matrix`, :func:`simple_cycles`,
  :func:`graph_relevant_vectors`, :func:`graph_face_counts`, :func:`tropical_moment`, :func:`graph_second_moment_function`

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
from VoronoiCellToolBox.graph_lattice import (
    graph_edges,
    cycle_basis_matrix,
    riemann_matrix,
    simple_cycles,
    graph_relevant_vectors,
    graph_face_counts,
    tropical_moment,
    graph_second_moment_function,
)

__version__ = "2.1.0"
__author__ = "Yelena Mandelshtam, Raul Penaguiao"

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
    "graph_edges",
    "cycle_basis_matrix",
    "riemann_matrix",
    "simple_cycles",
    "graph_relevant_vectors",
    "graph_face_counts",
    "tropical_moment",
    "graph_second_moment_function",
    # legacy names
    "VCell",
    "eq_of_wall",
    "Qform",
    "Delsets",
    "FormatPullingTrigMatrix",
    "smMatrixSimplex",
    "normalizedChamberSecondMomentPolynomial",
]
