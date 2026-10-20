r"""
Pictures of Voronoi cells in dimensions 2 and 3
"""
from sage.all import Graphics, polygon, line, point
from VoronoiCellToolBox.voronoi_cell import voronoi_cell, pulling_triangulation, _rows


def plot_voronoi_cell(Q, triangulation=False, relevant=False, **options):
    r"""
    Plot the Voronoi cell of `Q` (dimension 2 or 3).

    INPUT:

    - ``Q`` -- a positive definite metric of size 2 or 3
    - ``triangulation`` -- (default: ``False``) draw the pulling triangulation used by
      :func:`second_moment`; in dimension 3 its edges are drawn
    - ``relevant`` -- (default: ``False``) in dimension 2, mark the relevant vectors
    - further keyword arguments are passed to Sage's ``Polyhedron.plot``

    The cell is drawn in the coordinates of `\ZZ^d`, in which it is the region closer to `0`
    than to any other lattice point in the metric `Q` (so it looks sheared for non-diagonal `Q`).

    EXAMPLES::

        sage: from VoronoiCellToolBox import plot_voronoi_cell
        sage: g = plot_voronoi_cell([[2, -1], [-1, 2]], triangulation=True, relevant=True)
        sage: type(g)
        <class 'sage.plot.graphics.Graphics'>
        sage: g3 = plot_voronoi_cell([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]], triangulation=True)
    """
    Q = _rows(Q)
    d = len(Q)
    if d not in (2, 3):
        raise ValueError("only dimensions 2 and 3 can be plotted")
    V = voronoi_cell(Q)
    if d == 2:
        g = V.plot(fill="lightblue", alpha=0.6, **options)
        if triangulation:
            for simplex in pulling_triangulation(V):
                g += polygon([list(v) for v in simplex], fill=False, color="black", thickness=0.7)
        if relevant:
            from VoronoiCellToolBox.voronoi_cell import relevant_vectors
            g += point(relevant_vectors(Q), color="red", size=30)
            g += point([(0, 0)], color="black", size=30)
        g.set_aspect_ratio(1)
        return g
    g = V.plot(fill="lightblue", alpha=0.3, **options)
    if triangulation:
        g3 = Graphics()
        for simplex in pulling_triangulation(V):
            pts = [list(v) for v in simplex]
            for i in range(len(pts)):
                for j in range(i + 1, len(pts)):
                    g3 += line([pts[i], pts[j]], color="black", thickness=0.5)
        g += g3
    return g
