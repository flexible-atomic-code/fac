from distutils.core import setup, Extension
import os

# change the fortran runtime lib. appropriately.
# this only works for sun solaris.
fortranlib = ["F77", "M77", "sunmath"]

# for Linux, you would use:
# fortranlib = ["g77"]

libs = ["fac", "lapack", "blas", "coul",
        "toms", "mpfun", "minpack", "m"] + fortranlib

libdir = ["faclib", "lapack", "blas", "coul",
          "toms", "mpfun", "minpack"]
incdir = ["faclib"]
    
os.system("make lib")
os.system("make doc")

setup(name = "FAC",
      version = "0.6.7",
      package_dir = {'pfac': 'python'},
      py_modules = ['pfac.const', 'pfac.config', 'pfac.dr'],
      ext_modules = [Extension("pfac.fac",
                               ["python/fac.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs)])

