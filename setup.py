#
#   FAC - Flexible Atomic Code
#   Copyright (C) 2001-2015 Ming Feng Gu
# 
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#

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
            bsfac = 'SFAC-%s.%s.tar'%(version, get_platform())
      elif (s[0:4] == '-mpy'):
            pyinc = get_config_var('INCLUDEPY')
            pylib = get_config_var('LIBPL')+'/'+get_config_var('LDLIBRARY')
            os.system('cd python; export PYINC=%s PYLIB=%s; make mpy'%(pyinc, pylib))
            no_setup = 1
      else:
            sys.argv.append(s)

if (no_setup == 0):
      setup(name = "PFAC",
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
      print('Creating SFAC binary ...')
      c = 'cd sfac; mkdir bin;'
      c += 'cp sfac bin; cp scrm bin; cp spol bin;'
      c += 'tar cvf %s ./bin;'%bsfac
      c += 'gzip %s;'%bsfac
      c += 'rm -rf bin;'
      c += 'mv %s.gz ../dist/;'%bsfac
      os.system(c)
      
