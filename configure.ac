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
# actually any crtn.o may cause trouble, we don't need that.
AC_MSG_CHECKING([for -lcrtn.o in FLIBS])
FLIBS=[`echo $FLIBS | sed 's/-lcrt.*\.o//g'`]
FLIBS=[`echo $FLIBS | sed 's/  */ /g'`]
AC_MSG_RESULT([$FLIBS])

#get optimization option for C and f77 compilers
AC_ARG_ENABLE(copt,
	[ --enable-copt=opt	optimization option for C compiler],
	[copt="$enableval"])
AC_ARG_ENABLE(fopt,
	[ --enable-fopt=opt	optimization option for F77 Compiler],
	[fopt="$enableval"])

# get position independent code option for C and f77 compilers
AC_ARG_ENABLE(cpic,
	[ --enable-cpic=opt	PIC option for C compiler],
	[cpic="$enableval"])
AC_ARG_ENABLE(fpic,
	[ --enable-fpic=opt	PIC option for F77 Compiler],
	[fpic="$enableval"])

# override the default CFLAGS and FFLAGS for GCC and G77
if test -n "$GCC"
then
  if test -z "$copt"
  then
    copt="-Ofast"
  fi
  if test -z "$cpic"
  then
    cpic="-fPIC"
  fi
fi
if test -n "$G77"
then
  if test -z "$fopt"
  then
    fopt="-O2"
  fi
  if test -z "$fpic"
  then
    fpic="-fPIC"
  fi
fi
CFLAGS="$copt $cpic"
FFLAGS="$fpic"

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
  if test -n "$gcc_warnings"
  then
    CFLAGS="$CFLAGS -Wall -W -pedantic -Winline -Wmissing-prototypes \
 -Wnested-externs -Wpointer-arith -Wcast-align -Wshadow -Wstrict-prototypes"
    CFLAGS=`echo $CFLAGS`
  fi
fi
if test -n "$G77"
then
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
      *hpux*)
	xxFortran=hpuxFortran
      ;;
      *)
        echo "NOTE: F77 known by cfortran.h is not detected."
        echo "      using default f2cFortran."
        echo "      if not, change it using --with-cfortran=xxFortran"
        echo "      see faclib/cfortran.doc for available options."
        xxFortran=f2cFortran
      ;;
    esac
  fi
fi
AC_DEFINE_UNQUOTED([$xxFortran])

# MPI
AC_ARG_WITH(mpi,
	[ --with-mpi		Use MPI],
	[use_mpi=$withval])
AC_ARG_WITH(mpicompile,
	[ --with-mpicompile	MPI compile flags],
	[mpicompile=$withval])
AC_ARG_WITH(mpilink,
	[ --with-mpilink	MPI link flags],
	[mpilink=$withval])
AC_ARG_WITH(mpifflag,
	[ --with-mpifflag	MPI f77 flags],
	[mpifflag=$withval])
AC_ARG_WITH(extrainc,
	[ --with-extrainc	Extra compile flags],
	[extrainc=$withval])
AC_ARG_WITH(extralib,
	[ --with-extralib	Extra Libs],
	[extralib=$withval])
AC_ARG_WITH(blaslib,
	[ --with-blaslib	BLAS Libs],
	[blaslib=$withval])
AC_ARG_WITH(lapacklib,
	[ --with-lapacklib	LAPACK Libs],
	[lapacklib=$withval])

# setup libs
if test "x$blaslib" == "x"
then
  blaslib="-lmblas"
fi
if test "x$lapacklib" == "x"
then
  lapacklib="-lmlapack"
fi
LIBS="$LIBS -lfac $lapacklib $blaslib -lpthread"

fwithm=`echo "$FLIBS " | grep " -lm"`
if test "x$fwithm" = "x"
then
  LIBS="$LIBS -lm"
fi

if test "x$use_mpi" != "x"
then
  case "$use_mpi" in
    *mpi|1)
      if test "x$mpicompile" = "x"
      then
        mpicompile=`mpicc -showme`
        mpicompile=[`echo $mpicompile | sed 's/^ *//'`]
	mpicompile=[`echo $mpicompile | sed 's/^[^ ]* *//'`]
	mpicompile=[`echo $mpicompile | sed 's/ -[Ll][^ ]*//g'`]
	mpicompile=[`echo $mpicompile | sed 's/ *$//'`]
      fi
      if test "x$mpilink" = "x"
      then
        mpilink=`mpicc -showme`
	mpilink=[`echo $mpilink | sed 's/^ *//'`]
	mpilink=[`echo $mpilink | sed 's/^[^ ]* *//'`]
	mpilink=[`echo $mpilink | sed 's/ -[^Ll][^ ]*//g'`]
	mpilink=[`echo $mpilink | sed 's/ *$//'`]
      fi
      AC_DEFINE_UNQUOTED([USE_MPI])
    ;;
    omp|2)
      if test "x$mpicompile" = "x"
      then
  	mpicompile=-fopenmp
	mpilink="-fopenmp -lgcc_eh"
      fi
      if test "x$mpifflag" = "x"
      then
        mpifflag=-fopenmp
      fi
      AC_DEFINE_UNQUOTED([USE_MPI], [2])
    ;;
    *)
      AC_DEFINE_UNQUOTED([USE_MPI])
    ;;
  esac

  CPPFLAGS="$CPPFLAGS $mpicompile"
  LIBS="$LIBS $mpilink"
  FFLAGS="$FFLAGS $mpifflag"
fi

if test "x$extrainc" != "x"
then
  CPPFLAGS="$CPPFLAGS $extrainc"
fi
if test "x$extralib" != "x"
then
  LDFLAGS="$LDFLAGS $extralib"
fi

# special flags for some OS
case "$host_os" in
  *darwin*)
    if test -n "$GCC"
    then
      AC_CHECK_LIB(cc_dynamic,main,LIBS="$LIBS -lcc_dynamic")
    fi
  ;;
  *)
  ;;
esac

# add fopt
FFLAGSNO="$FFLAGS"
AC_SUBST(FFLAGSNO)
FFLAGS="$FFLAGS $fopt"

# add include pathes.
CPPFLAGS="$CPPFLAGS -I$TOPDIR"
CPPFLAGS="$CPPFLAGS -I$TOPDIR/faclib"
LDFLAGS="$LDFLAGS -L$TOPDIR"

AC_CONFIG_FILES([Makefile
		 sfac/Makefile
		 python/Makefile
		 faclib/Makefile
		 blas/Makefile
		 coul/Makefile
		 ionis/Makefile
		 lapack/Makefile
		 minpack/Makefile
		 mpfun/Makefile
		 ode/Makefile
		 toms/Makefile
		 modqed/Makefile
		 quadpack/Makefile])

if test -z "$PYTHON"
then
  echo "WARNING: NO PYTHON DETECTED ON THIS SYSTEM."
  echo "To install PFAC interface,"
  echo "download and install Python 2.0 or later"
  echo "You may still procede to install SFAC interface, though."
fi

AC_OUTPUT
