# Process this file with autoconf to produce a configure script.
AC_INIT

AC_CONFIG_HEADERS([sysdef.h:sysdef.h.in])
AC_PREFIX_DEFAULT("/usr/local")

TOPDIR=`pwd`
AC_SUBST(TOPDIR)

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_CHECK_PROG([PYTHON], [python], [python])
AC_SUBST(PYTHON)

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_HEADER_STDC

AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS

# some F77 LDFLAGS have a -lcrt1.o which causes trouble.
# actuall any crtn.o may cause trouble, we don't need that.
AC_MSG_CHECKING([for -lcrtn.o in FLIBS])
NFLIBS=[`echo $FLIBS | sed 's/-lcrt[0-9]\.o//g'`]
NFLIBS=[`echo $NFLIBS | sed 's/  */ /g'`]
if test "x$NFLIBS" = "x$FLIBS"
then
  AC_MSG_RESULT([no.])
else
  FLIBS=$NFLIBS
  AC_MSG_RESULT([$FLIBS])
fi

# override the default CFLAGS and FFLAGS
if test -n "$GCC"
then
  CFLAGS="-O2"
else
  CFLAGS="-O"
fi
if test -n "$G77"
then 
  FFLAGS="-O2"
else
  FFLAGS="-O"
fi

# Debugging option
AC_ARG_ENABLE(cdebug,
	[  --enable-cdebug=opt    Debugging option for C],
	[cdebug_opt=$enableval])
if test -n "$cdebug_opt"
then
  CFLAGS="$CFLAGS $cdebug_opt"
fi
AC_ARG_ENABLE(fdebug,
	[  --enable-fdebug=opt    Debugging option for F77],
	[fdebug_opt=$enableval])
if test -n "$fdebug_opt"
then
  FFLAGS="$FFLAGS $fdebug_opt"
fi
  
# gcc warnings
AC_ARG_ENABLE(warnings, 
	[  --enable-warnings	GCC warnings],
	[gcc_warnings=$enableval])
if test -n "$GCC"
then
  CFLAGS="$CFLAGS -fno-strength-reduce"
  if test -n "$gcc_warnings"
  then
    CFLAGS="$CFLAGS -Wall -W -pedantic -Winline -Wmissing-prototypes \
 -Wnested-externs -Wpointer-arith -Wcast-align -Wshadow -Wstrict-prototypes"
    CFLAGS=`echo $CFLAGS`
  fi
fi
if test -n "$G77"
then
  FFLAGS="$FFLAGS -fno-strength-reduce"
  if test -n "$gcc_warnings"
  then
    FFLAGS="$FFLAGS -Wall"
    FFLAGS=`echo $FFLAGS`
  fi
fi

# which F77 known by cfortran.h is used
AC_ARG_WITH(cfortran, 
	[  --with-cfortran=xxFortran    Macro to be defined for cfortran.h],
	[xxFortran="$withval"])
if test -z "$xxFortran"
then
  if test -n "$G77"
  then 
    xxFortran=f2cFortran
  else
    case "$host_os" in
      *linux* )
        xxFortran=f2cFortran
      ;;
      *solaris*)
	xxFortran=sunFortran
      ;;
      *sunos*)
	xxFortran=sunFortran
      ;;
      *)
        echo "NOTE: F77 known by cfortran.h is not detected."
        echo "      using default f2cFortran."
        echo "      if not, chage it using --with-cfortran=xxFortran"
        echo "      see faclib/cfortran.doc for available options."
        xxFortran=f2cFortran
      ;;
    esac
  fi
fi
AC_DEFINE_UNQUOTED([$xxFortran])

# special flags for some OS
case "$host_os" in
  *darwin*)
    CFLAGS="$CFLAGS -no-cpp-precomp -fno-common"
    if test -n "$G77"
    then 
      FFLAGS="$FFLAGS -fno-common"
    fi
  ;;
  *)
  ;;
esac 
      
# add include and lib. pathes.	      
CPPFLAGS="$CPPFLAGS -I$TOPDIR"
CPPFLAGS="$CPPFLAGS -I$TOPDIR/faclib"
LDFLAGS="$LDFLAGS -L$TOPDIR"
LIBS="$LIBS -lfac"
fwithm=`echo "$FLIBS " | grep " -lm"`
if test "x$fwithm" = "x"
then
  LIBS="$LIBS -lm"
fi

AC_CONFIG_FILES([Makefile 
		 sfac/Makefile
		 faclib/Makefile 
		 blas/Makefile 
		 coul/Makefile
		 ionis/Makefile
		 lapack/Makefile
		 minpack/Makefile 
		 mpfun/Makefile
		 ode/Makefile
		 toms/Makefile
		 quadpack/Makefile])

if test -z "$PYTHON"
then
  echo "WARNING: NO PYTHON DETECTED ON THIS SYSTEM."
  echo "To install PFAC interface,"
  echo "download and install Python 2.0 or later"
  echo "You may still procede to install SFAC interface, though."
fi

AC_OUTPUT