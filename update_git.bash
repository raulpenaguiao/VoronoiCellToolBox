git add setup.py update_git.bash VoronoiCellToolBox/orchestration.py VoronoiCellToolBox/templatecomputation.m2 demo.ipynb VoronoiCellToolBox/macaulay_parsing.py VoronoiCellToolBox/voronoi_cell.py
git commit -m "Updated to version 1.3.6"
git push
sage --pip install --force git+https://github.com/raulpenaguiao/VoronoiCellToolBox
