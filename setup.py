"""
Setup configuration for MolOrbDraw package.
"""
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="molorbdraw",
    version="0.1.0",
    author="MolOrbDraw Contributors",
    description="A tool for visualizing molecular orbitals from cube files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PhelanShao/molorbdraw",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "vtk>=9.0.0",
    ],
    entry_points={
        "console_scripts": [
            "molorbdraw=molorbdraw.main:main",
        ],
    },
)
