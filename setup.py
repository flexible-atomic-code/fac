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
from distutils import sysconfig
import distutils
from distutils.ccompiler import CCompiler
from distutils.unixccompiler import UnixCCompiler
from distutils.util import get_platform
import sys
import os

incdir = []
libdir = []
libs = []
extralink = []
extracomp = []
bsfac = ''

no_setup = 0
x = sys.argv[2:]
sys.argv = sys.argv[0:2]

CC = None


# Obtain version & consts from faclib/consts.h
def parse_consts(filename=None, vfile=None, cfile=None):
    if not filename:
        filename = os.path.join(
            os.path.dirname(__file__), 'faclib', 'consts.h')
    if not vfile:
        vfile = os.path.join(
            os.path.dirname(__file__), 'python', 'version.py')
    if not cfile:
        cfile = os.path.join(
            os.path.dirname(__file__), 'python', 'consts.py')
    consts = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            a = line.strip().split()
            if len(a) < 3:
                continue
            if a[0] != '#define':
                continue
            consts[a[1]] = a[2]
    a = open(vfile, 'w')
    version = '%s.%s.%s'%(consts['VERSION'], consts['SUBVERSION'], consts['SUBSUBVERSION'])
    try:
        a.write("version = '%s'\n"%version)
    finally:
        a.close()
    a = open(cfile, 'w')
    try:
        for k in consts.keys():
            v = consts[k]
            try:
                if v[0] == '(' and v[-1] == ')':
                    v = v[1:-1]
                dv = float(v)
                a.write("%s = %s\n"%(k, v))
            except:
                pass
    finally:
        a.close()
    return version

VERSION = parse_consts()

for s in x:
    if (s[0:11] == '-ccompiler='):
        CC = s[11:]
    elif (s[0:11] == '-extracomp='):
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
        bsfac = 'SFAC-%s.%s.tar'%(VERSION, get_platform())
    elif (s[0:4] == '-mpy'):
        pyinc = sysconfig.et_config_var('INCLUDEPY')
        pylib = (sysconfig.get_config_var('LIBPL') +
                 '/' + sysconfig.get_config_var('LDLIBRARY'))
        os.system('cd python; export PYINC=%s PYLIB=%s; make mpy'%(pyinc, pylib))
        no_setup = 1
    else:
        sys.argv.append(s)


# We studied a lot from Intel/pyMIC repoif CC is None:
# https://github.com/intel/pyMIC/blob/master/setup.py
# compiler driver for Anaconda Composer for C/C++
class MyCompiler(UnixCCompiler, object):
    """Compiler wrapper for anaconda gcc_linux-64 """
    def set_executables(self, **args):
        # basically, we ignore all the tool chain coming in
        if CC is not None:
            super(self.__class__, self).set_executables(
                compiler=CC, compiler_so=CC, linker_exe=CC,
                linker_so=CC)

    def _fix_lib_args(self, libraries, library_dirs, runtime_library_dirs):
        # we need to have this method here, to avoid an endless
        # recursion in UnixCCompiler._fix_lib_args.
        libraries, library_dirs, runtime_library_dirs = \
            CCompiler._fix_lib_args(self, libraries, library_dirs,
                                    runtime_library_dirs)
        libdir = sysconfig.get_config_var('LIBDIR')
        if runtime_library_dirs and (libdir in runtime_library_dirs):
            runtime_library_dirs.remove(libdir)
        return libraries, library_dirs, runtime_library_dirs


# register Intel Composer for C/C++ as the main compiler for pyMIC
def my_new_compiler(plat=None, compiler=None, verbose=0, dry_run=0, force=0):
    compiler = MyCompiler(verbose, dry_run, force)
    return compiler


# override compiler construction
distutils.ccompiler.new_compiler = my_new_compiler

if CC is not None:
    ldshared = sysconfig.get_config_var('LDSHARED')
    ldshared = ldshared.split(' ')
    for a in ldshared[1:]:
        extralink.append(a)
    
Extensions = [Extension(mod, src, include_dirs=incdir, library_dirs=libdir,
                        libraries=libs, extra_compile_args=extracomp,
                        extra_link_args=extralink)
              for mod, src in (("pfac.fac", ["python/fac.c"]),
                               ("pfac.crm", ["python/pcrm.c"]),
                               ("pfac.pol", ["python/ppol.c"]),
                               ("pfac.util", ["python/util.c"]))]

if (no_setup == 0):
    setup(name="PFAC",
          version=VERSION,
          package_dir={'pfac': 'python'},
          py_modules=['pfac.const', 'pfac.config', 'pfac.table',
                      'pfac.atom', 'pfac.spm', 'pfac.rfac',
                      'pfac.version', 'pfac.consts'],
          ext_modules=Extensions)

if (sys.argv[1][0:5] == 'bdist' and bsfac != ''):
    print('Creating SFAC binary ...')
    c = 'cd sfac; mkdir bin;'
    c += 'cp sfac bin; cp scrm bin; cp spol bin;'
    c += 'tar cvf %s ./bin;'%(bsfac)
    c += 'gzip %s;'%(bsfac)
    c += 'rm -rf bin;'
    c += 'mv %s.gz ../dist/;'%(bsfac)
    os.system(c)
