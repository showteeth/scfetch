# Installation

Here, we list some issues related to installation.

## C++14/C++17 standard requested but CXX14/CXX17 is not defined

This is a very common issue about **gcc version** (the following is solution on `CentOS 7`):

* Install gcc (root required):
```bash
# check gcc
gcc --version
  # gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-44)

# install gcc
yum install centos-release-scl
yum install devtoolset-8
scl enable devtoolset-8 bash

# check gcc version
gcc --version
  # gcc (GCC) 8.3.1 20190311 (Red Hat 8.3.1-3)
```

* Configure R env:
```bash
# create ~/.R/Makevars first
# change CXX to installed gcc path
cat ~/.R/Makevars
  # MAKEFLAGS = -j8
  # 
  # ## C++ flags
  # CXX=/opt/rh/devtoolset-8/root/usr/bin/g++
  # CXX11=/opt/rh/devtoolset-8/root/usr/bin/g++
  # CXX14=/opt/rh/devtoolset-8/root/usr/bin/g++
  # CXX17=/opt/rh/devtoolset-8/root/usr/bin/g++
  # 
  # CXXFLAGS=-O3 -march=native -Wno-ignored-attributes
  # CXX11FLAGS=-O3 -march=native -Wno-ignored-attributes
  # CXX14FLAGS=-O3 -march=native -Wno-ignored-attributes
  # CXX17FLAGS=-O3 -march=native -Wno-ignored-attributes
  # 
  # CXXPICFLAGS=-fPIC
  # CXX11PICFLAGS=-fPIC
  # CXX14PICFLAGS=-fPIC
  # CXX17PICFLAGS=-fPIC
  # 
  # CXX11STD=-std=c++11
  # CXX14STD=-std=c++14
  # CXX17STD=-std=c++17
  # 
  # ## C flags
  # CC=gcc
  # CFLAGS=-O3 -march=native
  # 
  # ## Fortran flags
  # FC=gfortran
  # F77=gfortran
  # FFLAGS=-O3 -march=native
  # FCFLAGS=-O3 -march=native
```

[Reference](https://github.com/stan-dev/rstan/issues/892)

<hr />

## hdf5: Please ensure that at least version 1.8.13 is installed

CentOS 7 installs hdf5 `1.8.12` by default. The following is solution:

```bash
# download hdf5: https://www.hdfgroup.org/downloads/hdf5/source-code/

# cofigure
tar -zxvf hdf5-1.14.2.tar.gz
cd hdf5-1.14.2/
./configure --prefix=/usr/local/hdf5
make
make install
```

Install:

```R
install.packages("hdf5r", configure.args="--with-hdf5=/usr/local/hdf5/bin/h5cc")
```

[Reference](https://github.com/hhoeflin/hdf5r/issues/115)

<hr />

## SeuratObject

```R
devtools::install_version("SeuratObject", version = "4.1.4", repos = "https://cran.r-project.org")
```

<hr />

## General solution

For general purpose, you can create a new R environment (same R version) using [Conda](), install packages via Conda, add R library path of Conda to `.libPath()`.

* Install Conda: [reference](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
* Install R: `conda install r-base=4.0.3` (keep the same version)
* Install packages: `conda install xxxx`
* Add to R path: `.libPaths(c(.libPaths(), "/home/user/R/x86_64-conda-linux-gnu-library/4.0"))`

This can solve most of the problems!!

<hr />



