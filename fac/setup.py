from distutils.core import setup, Extension
import os
import sys

# set the fortran runtime lib. appropriately.
# this only works for sun solaris.
# fortranlib = ["F77", "M77", "sunmath"]
# for Linux (or Windows with Cygwin), you would use:
# fortranlib = ["g2c"]

platform = sys.platform[:5]
if (platform == 'linux' or
    platform == 'cygwi'):
      fortranlib = ["g2c"]
      fflags = '-O -c -fPIC'
      cflags = '-O -c -fPIC'
elif (platform == 'sunos'):
      fortranlib = ["F77", "M77", "sunmath"]
      fflags = '-O -c -KPIC'
      cflags = '-O -c -KPIC'

flib = ' '.join(map(lambda a:'-l%s'%a, fortranlib))
macros = '"FORTRANLIB=%s" '%flib
macros = macros + '"FFLAGS=%s" '%fflags
macros = macros + '"CFLAGS=%s" '%cflags
print macros
os.system("make lib %s"%macros)
os.system("make doc")

libs = ["fac", "quadpack", "lapack", "blas", "coul",
        "toms", "mpfun", "minpack", "ionis", "m"] + fortranlib
libdir = ["faclib", "coul", "toms", "mpfun", "minpack",
          "quadpack", "lapack", "blas", "ionis"]
incdir = ["faclib"]
    
setup(name = "FAC",
      version = "0.7.1",
      package_dir = {'pfac': 'python'},
      py_modules = ['pfac.const', 'pfac.config', 'pfac.atom', 'pfac.spm'],
      ext_modules = [Extension("pfac.fac",
                               ["python/fac.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs),
                     Extension("pfac.crm",
                               ["python/pcrm.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs)
                     ])
