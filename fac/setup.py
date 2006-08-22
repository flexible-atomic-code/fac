from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var
from distutils.util import get_platform
import sys
import os

incdir = []
libdir = []
libs = []
extralink = []
extracomp = []
bsfac = ''
version = "1.1.1"

no_setup = 0
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
      elif (s[0:6] == '-bsfac'):
            bsfac = 'SFAC-%s.%s.tar.gz'%(version, get_platform())
      elif (s[0:4] == '-mpy'):
            pyinc = get_config_var('INCLUDEPY')
            pylib = get_config_var('LIBPL')+'/'+get_config_var('LDLIBRARY')
            os.system('cd python; export PYINC=%s PYLIB=%s; make mpy'%(pyinc, pylib))
            no_setup = 1
      else:
            sys.argv.append(s)

if (no_setup == 0):
      setup(name = "FAC",
            version = version,
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

if (sys.argv[1][0:5] == 'bdist' and bsfac != ''):
      print 'Creating SFAC binary ...'
      os.system('cd sfac; tar zcvf %s sfac scrm spol'%bsfac)
      os.system('mv sfac/%s dist'%bsfac)
      
