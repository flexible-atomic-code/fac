from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var
import sys
import os

incdir = []
libdir = []
libs = []
extralink = []
extracomp = []

x = sys.argv[2:]
sys.argv = sys.argv[0:2]
for s in x:
      if (s[0:11] == '-extracomp='):
            for a in s[11:].split():
                  if (a[0:2] == '-I'):
                        incdir.append(a[2:])
                  elif (a[0:2] == '-L'):
                        libdir.append(a[2:])
                  elif (a[0:2] == '-l'):
                        libs.append(a[2:])
                  else:
                        extracomp.append(a)
      elif (s[0:11] == '-extralink='):
            for a in s[11:].split():
                  if (a[0:2] == '-I'):
                        incdir.append(a[2:])
                  elif (a[0:2] == '-L'):
                        libdir.append(a[2:])
                  elif (a[0:2] == '-l'):
                        libs.append(a[2:])
                  else:
                        extralink.append(a)
      else:
            sys.argv.append(s)

setup(name = "FAC",
      package_dir = {'pfac': 'python'},
      py_modules = ['pfac.const', 'pfac.config', 'pfac.table',
                    'pfac.atom', 'pfac.spm'],
      ext_modules = [Extension("pfac.fac",
                               ["python/fac.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs,
                               extra_compile_args = extracomp,
                               extra_link_args = extralink),
                     Extension("pfac.crm",
                               ["python/pcrm.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs,
                               extra_compile_args = extracomp,
                               extra_link_args = extralink),
                     Extension("pfac.pol",
                               ["python/ppol.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs,
                               extra_compile_args = extracomp,
                               extra_link_args = extralink),
                     Extension("pfac.util",
                               ["python/util.c"],
                               include_dirs = incdir,
                               library_dirs = libdir,
                               libraries = libs,
                               extra_compile_args = extracomp,
                               extra_link_args = extralink)
                     ]) 

if (sys.argv[1] == 'build'):
      pyinc = get_config_var('INCLUDEPY')
      pylib = get_config_var('LIBPL')+'/'+get_config_var('LIBRARY')
      os.system('cd python; export PYINC=%s PYLIB=%s; make mpy'%(pyinc, pylib))
      
