from distutils.core import setup, Extension
import sys

libs = []
libdir = ["."]
incdir = [".", "faclib"]
extralink = []

flibs = sys.argv[2:]
sys.argv = sys.argv[0:2]
for s in flibs:
      if (s[0:2] == '--'):
            sys.argv.append(s)
      else:
            extralink.append(s)

setup(name = "FAC",
      package_dir = {'pfac': 'python'},
      py_modules = ['pfac.const', 'pfac.config', 'pfac.table',
                    'pfac.atom', 'pfac.spm'],
      ext_modules = [Extension("pfac.fac",
                               ["python/fac.c"],
                               include_dirs = incdir, 
                               library_dirs = libdir,
                               extra_link_args = extralink,
                               libraries = libs),
                     Extension("pfac.crm",
                               ["python/pcrm.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               extra_link_args = extralink,
                               libraries = libs),
                     Extension("pfac.util",
                               ["python/util.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               extra_link_args = extralink,
                               libraries = libs)
                     ]) 

