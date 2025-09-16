import setuptools
import numpy as np
import itertools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='VoronoiCellToolBox',
    version="1.2.14",
    author="Yelena Mandelshtam, Raul Penaguiao",
    author_email="raulpenaguiao@proton.me",
    license="GPL2+",
    description="Evaluate Voronoi Cells in Sagemath using python and macaulay2",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    package_data={
        'VoronoiCellToolBox': ['templatecomputation.m2'],  # Include the template
    },
    include_package_data=True,
    include_dirs=[np.get_include()],
    zip_safe=False,
)