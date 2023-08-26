[![Travis CI](https://travis-ci.org/WillemEerland/teetool.svg?branch=master)](https://travis-ci.org/WillemEerland/teetool)

A modification version of teetool for python3.9 and for mat data analysis and plotting
updated in 20230826 by zhengcong

original version readme
# teetool
a package to support with the statistical analysis of trajectory data -- helpful at determining the probability of clusters (collection) of trajectories colliding

purely spatial, ignores temporal effects

publication is available at http://doi.org/10.5334/jors.163

documentation is available at https://willemeerland.github.io/teetool/

# setup the environment in Linux

- download & install Anaconda from https://www.continuum.io/download
- open terminal
- navigate to Teetool directory

> conda create -n teetool python=2.7 pytest pytest-cov mayavi numpy scipy matplotlib pyside

> source activate teetool

> pip install .

# setup the environment in macOS

- download & install Anaconda from https://www.continuum.io/download
- open terminal
- navigate to Teetool directory

> conda create -n teetool python=2.7 pytest pytest-cov mayavi numpy scipy matplotlib

> source activate teetool

> pip install .

# setup the environment in Windows

- download & install Anaconda from https://www.continuum.io/download
- open 'Anaconda prompt'
- navigate to Teetool directory

> conda create -n teetool python=2.7 pytest pytest-cov mayavi numpy scipy matplotlib pyside

> activate teetool

> set QT_API=pyside

> pip install .

# run tests

> py.test

# run tests, including coverage report

> (cd test ; py.test -v --cov-report html --cov=teetool)

# run examples via Jupyter notebook

> pip install .[example]

> jupyter notebook

- find example/ in browser and run files in order
