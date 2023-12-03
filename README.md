# PyHub
----------

PyHub is an easy-to-use, hybrid Python/Fortran package to solve efficiently the Fermi-Hubbard Model

Installation
------------

To install, clone the repository

```
git clone https://github.com/Quentin-Marecat/PyHub.git
```

Install HDF5-fortran BLAS and LAPACK libraries

```
sudo apt install  libhdf5-dev libhdf5-fortran-102 libopenblas-dev liblapack-dev
```

Install the package using `pip` from the top-level directory, which requires CMake

```
python -m pip install . --user
```

Add `pyhub` in your .bashrc

```
export PYTHONPATH="${PYTHONPATH}:/path/you/install/PyHub/"
```

Quickstart
----------

Examples of how to use PyHub can be found in the `pyhub/examples` directory.

Authors
----------

Q. Marécat
