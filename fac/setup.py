from distutils.core import setup, Extension
import os
import sys

# set the fortran runtime lib. appropriately.
# and the compiler flags. 

platform = sys.platform[:5]
if (platform == 'linux'):
      fortranlib = ['g2c']
      fflags = '-O -c -fPIC'
      cflags = '-O -c -fPIC'
      ldflags = ''
      cc = 'gcc'
      fc = 'g77'
elif (platform == 'cygwi'):
      fortranlib = ['g2c']
      fflags = '-O -c'
      cflags = '-O -c'
      ldflags = ''
      cc = 'gcc'
      fc = 'g77'
elif (platform == 'sunos'):
      fortranlib = ['F77', 'M77', 'sunmath']
      fflags = '-O -c -KPIC'
      cflags = '-O -c -KPIC'
      ldflags = ''
      cc = 'cc'
      fc = 'f77'
elif (platform == 'darwi'):
      fortranlib = ['g2c']
      fflags = '-O -c -fno-common'
      cflags = '-O -c -no-cpp-precomp -fno-common'
      ldflags = ''
      cc = 'cc'
      fc = 'f77'

flib = ' '.join(map(lambda a:'-l%s'%a, fortranlib))
macros = '"CC=%s" "FC=%s" '%(cc, fc)
macros = macros + '"FORTRANLIB=%s" '%flib
macros = macros + '"CFLAGS=%s" '%cflags
macros = macros + '"FFLAGS=%s" '%fflags 
macros = macros + '"LDFLAGS=%s" '%ldflags 

os.system("make lib %s"%macros)
if (sys.argv[1] == 'sfac'):
      os.system("make sfac %s"%macros)
      os.system("make install")
else:
      libs = ["fac", "fquadpack", "fode", "flapack", "fblas", "coul",
              "toms", "fmpfun", "fminpack", "ionis", "m"] + fortranlib
      libdir = ["faclib", "coul", "ode", "toms", "mpfun", "minpack",
                "quadpack", "lapack", "blas", "ionis"]
      incdir = ["faclib"]
    
      setup(name = "FAC",
            version = "0.8.4",
            package_dir = {'pfac': 'python'},
            py_modules = ['pfac.const', 'pfac.config', 'pfac.table',
                          'pfac.atom', 'pfac.spm'],
            ext_modules = [Extension("pfac.fac",
                                     ["python/fac.c"],
                                     include_dirs = incdir,
                                     library_dirs = libdir,
                                     libraries = libs),
                           Extension("pfac.crm",
                                     ["python/pcrm.c"],
                                     include_dirs = incdir,
                                     library_dirs = libdir,
                                     libraries = libs),
                           Extension("pfac.util",
                                     ["python/util.c"],
                                     include_dirs = incdir,
                                     library_dirs = libdir,
                                     libraries = libs)
                           ])
