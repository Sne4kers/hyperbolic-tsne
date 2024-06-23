# Accelerating hyperbolic t-SNE

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the code for the paper on application of the Cartesian quadtree for Barnes-Hut approximation for hyperbolic t-SNE. This is a fork of the original repository on the polar quadtree implementation distributed under MIT license and contains code unrelated to the Cartesian quadtree.

## Setup

Perform the following steps:

1. Install conda (we recommend using [miniconda](https://docs.conda.io/projects/miniconda/en/latest/))
2. Create environment: `conda create --name=htsne python=3.9.16`
3. Activate environment: `conda activate htsne`
4. Install dependencies with pip: `pip install -r requirements.txt`
5. Build Cython extensions: `python setup.py build_ext --inplace`
6. Install hyperbolic-tsne package: `pip install .`
7. To test installation run `python -c "from hyperbolicTSNE import HyperbolicTSNE"`. No errors should be raised and you should see the output `Please note that 'empty_sequence' uses the KL divergence with Barnes-Hut approximation (angle=0.5) by default.`.
8. To experiments and pictures from the paper, run scripts from `experiments_and_plots`.

Note 1: 
On macOS, the build process of the Cython extensions might yield an error if it cannot find OpenMP.
This error can be ignored and the package will still be correctly installed and able to run. 
The main consequence of this error is that the optimization iterations run slower.

## Use

In order to run either polar, cartesian `polar_or_cartesian="polar"` should be set to `"polar"` or `"cartesian"` respectively.
Look at the examples in `code.py`.

## Data

You can run hyperbolic TSNE on your high-dimensional data. 
Nevertheless, the examples and experiments in this repository rely on specific datasets. 
Below, we provide download links for each. 
We recommend putting all datasets in a `datasets` directory at the root of this repository.
The `load_data` function expects this path (`data_home`) to resolve the dataset.

Individual instructions per dataset:
- [LUKK](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-62)
- [MYELOID8000](https://github.com/scverse/scanpy_usage/tree/master/170430_krumsiek11)
- [PLANARIA](https://shiny.mdc-berlin.de/psca/)
- [MNIST](https://yann.lecun.com/exdb/mnist/)
- [WORDNET](https://github.com/facebookresearch/poincare-embeddings)
- [C_ELEGANS](https://github.com/Munfred/wormcells-data/releases)

## First steps

There are two ways of getting started with the `hyperbolicTSNE` package.
First, `code.py` offers a step-by-step guide showing how to use the HyperbolicTSNE package to embed a high-dimensional dataset.

## Replicating the paper results

This folder contains three types of files:
- Scripts to generate experimental data via embedding different data sets into hyperbolic space. These are pre-fixed with "data generation". 
- Scripts to create plots from the data, as they appear in the publication.
- Scripts to create tables from the data, as they appear in the publication.

The general workflow to reproduce the results from the paper is:
- Run the scripts to generate data.
- Run the scripts to plot the data.
- Run the scripts to generate tables.

Note that the data generation scripts assume a top-level folder, i.e., a folder next to "examples", "experiments", etc., called "datasets" that holds the datasets to be embedded.

## License and third-party software
The source code in this repository is released under the MIT License. However, all used third-party software libraries are governed by their respective licenses. Without the following libraries, this project would have been considerably harder: 
[scipy](https://scipy.org),
[numpy](https://numpy.org),
[scikit-learn](https://scikit-learn.org/stable/),
[hnswlib](https://github.com/nmslib/hnswlib),
[pandas](https://pandas.pydata.org),
[anndata](https://anndata.readthedocs.io/en/latest/),
[seaborn](https://seaborn.pydata.org),
[setuptools](https://github.com/pypa/setuptools),
[Cython](https://cython.org),
[tqdm](https://github.com/tqdm/tqdm),
[ipykernel](https://ipython.org).
