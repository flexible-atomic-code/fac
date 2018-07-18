# Flexible Atomic Code

[![travis](https://travis-ci.org/flexible-atomic-code/fac.svg?branch=master)](https://travis-ci.org/flexible-atomic-code/fac)

A software packge for the calculation of various atomic processes.

# Installation

There are two interfaces available to the FAC numerical lib. One makes use
of the Python scripting language, implemented as a package "pfac".
The other interface uses a simple command parser, compiled into executables.
The two interfaces are independent. Either one or both can be installed.


## 1. Download
Download the latest source packages from
[here](/https://github.com/flexible-atomic-code/fac/releases),
gunzip, and extract it.

## 2. Configure
Execute `configure` command.
You may need to specify some of the following options depending on your
environment.

#### 2-1. C/fortran compiler
If your C compiler is not gcc, then you need to supply the PIC (position
independent code) option with `--enable-cpic=***`.
If your F77 compiler is not g77, then you need to supply the PIC
option with `--enable-fpic=***`

Supply the optimization option for C and F77 compilers with
`--enable-copt=***` and `--enable-fopt=***` respectively if desired.

#### 2-2. cfortran
If your F77 compiler is not g77 or sun's f77, specify the option
`--with-cfortran=xxFortran`
see [faclib/cfortran.doc](faclib/cfortran.doc) for supported F77 compilers.

#### 2-3. Installation directory
If the default /usr/local is not what you want, specify `--prefix=my/dir`.
This dir only affects the installation of SFAC executables.

#### 2-4. Parallel computation
Some of the functions have a parallel version using MPI. You can build with MPI
enabled using the option `--with-mpi=***`, where `***` is the MPI implementation
installed on your machine. It has been tested with lammpi and openmp.
If you are using openmp, specify `--with-mpi=omp`.

If a different version of MPI is used, you have to supply the compile and link
flags to the C compiler with `--with-mpicompile` and `--with-mpilink`.


## 3. Install SFAC interface
To installs the SFAC interface, do
```
make
make install
```

## 4. Install Python interface
If you have Python 2.7 or later installed, then:
```
make pfac
make install-pfac
```
This installs the PFAC interface into Python's default site-package dir.

#### 4-1 Anaconda distribution
If your are using Anaconda distribution for Python environment, it may be
necessary to install gcc compilers from Anacona.
Linux:
+ gcc_linux-64
+ gfortran_linux-64  

macOS:
+ clang_osx-64
+ clangxx_osx-64
```
conda install gcc_linux-64 gfortran_linux-64
```
For the details of the Anaconda's compiler, see [the official page](https://conda.io/docs/user-guide/tasks/build-packages/compiler-tools.html).

It would be better to create a new virtual environment for FAC, as the compiler
installation sometimes affects other Python packages.
For the details of the virtual environment, see [Anaconda's documentation](https://conda.io/docs/user-guide/tasks/manage-python.html#installing-a-different-version-of-python)


# Usage

### PFAC interface

The modules can be called just like any python modules, either from the
interpreter interactively, or pass the script to python.
e.g.,
```
python fe24_level.py
```
Some example scripts can be found in [demo](demo) directory.

Note that PFAC is not thread safe. Only one instance of calculation can be done
in one Python kernel.


### SFAC Interface

The non-python interface produces 3 executables, "sfac", "scrm" and "spol",
accepting input files as arguments. The syntax of the input file for both
programs is identical, and is very similar to the Python syntax, except for
the missing of flow controls.

For more details, please read the [manual](doc/manual.pdf), and be sure to read
the FAQ section of the manual.


# Contents of FAC

### Third party content that cannot be distributed under GPL.

These directories are bundled in fac_util.tar.gz,
and are listed in the UtilList file.
+ blas/,       BLAS routines in fortran.
+ coul/,       fortran routines related to coulomb functions.
+ ionis/,      ionization balance calculation of Arnaud&Raymond,
and that of Mazzotta et al. Rates fitted by Verner et al.
+ lapack/,     fortran routines from LAPACK.
+ minpack/,    minpack optimization package.
+ quadpack/,   some routines from quadpack for quadrature.
+ mpfun/,      it used to be the multi-precision floating point arithmetic
package by Beiley. Since v0.8.6, it was replaced by my own implementation of
quad-precision arithmetics. The associated Legendre functions are also put into
this directory.
+ ode/,	     Livermore Solver for odinary differential equations.
+ toms/,       some numerical routines from TOMS.

## Contents distributed under GPL.

These files are bundled in fac_core.tar.gz
+ README,      this file.
+ ChangeLog,   changes.
+ configure.ac,AutoConf input file
+ configure,   configure script after parsing configure.ac
+ setup.py,    python setup script for building the pfac interface
+ demo/,       some demo scripts.
+ doc/,        reference guide for FAC
+ doc/papers/, papers describe the theoretical backgrounds.
+ faclib/,     core numerical library.
+ python/,     python interface.
+ sfac/,       non-python interface.


# Contribution

For bugs and suggestions, please report to
[our issue page](https://github.com/flexible-atomic-code/fac/issues).


# Contact

mfgu@ssl.berkeley.edu
Ming Feng Gu.
Space Science Laboratory
University of California Berkeley
