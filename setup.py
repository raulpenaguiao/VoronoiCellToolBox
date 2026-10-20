import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='VoronoiCellToolBox',
    version="2.1.0",
    author="Yelena Mandelshtam, Raul Penaguiao",
    author_email="raulpenaguiao@proton.me",
    license="GPL-3.0-or-later",
    description="Voronoi cells, secondary cones and exact chamber-wise second moments of lattices, in SageMath and Macaulay2",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/raulpenaguiao/VoronoiCellToolBox",
    packages=["VoronoiCellToolBox"],
    package_data={
        'VoronoiCellToolBox': ['*.m2'],  # the Macaulay2 code loaded at run time
    },
    include_package_data=True,
    python_requires=">=3.9",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    zip_safe=False,
)
