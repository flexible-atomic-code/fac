from distutils.core import setup, Extension
import sys

libs = ["fac", "m"]
libdir = ["."]
incdir = [".", "faclib"]
extralink = []

flibs = sys.argv[2:]
sys.argv = sys.argv[0:2]
nt = 'arg'
for s in flibs:
      s.strip()
      if (s[0:2] == '-L'):
	    if (len(s) > 2):
                  libdir.append(s[2:])
            else:
                  nt = 'dir'
      elif (s[0:2] == '-l'):
            if (len(s) > 2):
                  libs.append(s[2:])
            else:
                  nt = 'lib'
      else:
            if (nt == 'dir'):
                  libdir.append(s)
                  nt = 'arg'
            elif (nt == 'lib'):
                  libs.append(s)
            elif (nt == 'arg'):
                  if (s[0:1] == '-' and s[1:2] != '-'):
                        if (len(s) > 2):
                              extralink.append(s)
                        else:
                              nt = s[0:2]
                  elif (s[0:2] != '--'):
                        extralink.append(s)
                  else:
                        sys.argv.append(s)
            else:
                  extralink.append(nt+' '+s)
                  nt = 'arg'

setup(name = "FAC",
      version = "0.8.8",
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

