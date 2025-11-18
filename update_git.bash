git add setup.py update_git.bash 
git add VoronoiCellToolBox/orchestration.py VoronoiCellToolBox/templatecomputation.m2 
git add demo_secondmoment.ipynb demo_interfaceM2.ipynb
git add VoronoiCellToolBox/macaulay_parsing.py VoronoiCellToolBox/voronoi_cell.py
read -r -p "Enter commit message: " COMMIT_MSG
if [ -z "$COMMIT_MSG" ]; then
    echo "No commit message provided. Aborting."
    exit 1
fi
git commit -m "$COMMIT_MSG"
git push
sage --pip install --force git+https://github.com/raulpenaguiao/VoronoiCellToolBox
