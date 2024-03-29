# Process this file with autoconf to produce a configure script.
#
#
#
#

define([psi_version],[3.4.0])
define([psi_buildid],[alpha])
define([psi_bugreport],[psicode@users.sourceforge.net])
AC_PREREQ([2.71])
AC_INIT([psi],[psi_version],[psi_bugreport])

AC_CONFIG_HEADERS([include/psiconfig.h])

PSI_VERSION=psi_version
PSI_BUILDID=psi_buildid
PSI_BUGREPORT=psi_bugreport
AC_SUBST(PSI_VERSION)
AC_SUBST(PSI_BUILDID)
AC_SUBST(PSI_BUGREPORT)
AC_CONFIG_AUX_DIR(bin)

AC_CANONICAL_TARGET
AC_PREFIX_DEFAULT(/usr/local/psi)

AC_PROG_MAKE_SET()
# Psi now is primarily C++
AC_LANG([C++])

srcdir=`(cd $srcdir; pwd)`


objincdir='$(top_objdir)/include'
AC_SUBST([objincdir])

if test -z $docdir; then
       docdir='$(datadir)/doc/$(PACKAGE)'
       AC_SUBST([docdir])
fi
if test -z $txtdir; then
       txtdir='$(docdir)/txt'
       AC_SUBST([txtdir])
fi
if test -z $htmldir; then
       htmldir='$(docdir)/html'
       AC_SUBST([htmldir])
fi
if test -z $dvidir; then
       dvidir='$(docdir)/dvi'
       AC_SUBST([dvidir])
fi
pkgdatadir='$(datadir)/psi'
AC_SUBST([pkgdatadir])

INSTALL_INCLUDE='$(top_srcdir)/bin/installinc.sh' AC_SUBST(INSTALL_INCLUDE)
MKDIRS='$(top_srcdir)/bin/mkdirs.sh' AC_SUBST(MKDIRS)

########
## PSI-specific options
########

# Options for libint
LIBINT_NEW_AM=8
LIBINT_OPT_AM=4
AC_ARG_WITH(max-am-eri,[  --with-max-am-eri       maximum AM in ERIs],[
  LIBINT_NEW_AM=$((2*($withval)))
  LIBINT_OPT_AM=$(($withval))
])
AC_SUBST(LIBINT_NEW_AM)
AC_SUBST(LIBINT_OPT_AM)

# Options for libderiv
LIBDERIV_NEW_AM1=6
AC_ARG_WITH(max-am-deriv1,[  --with-max-am-deriv1    maximum AM in ERI first derivatives],[
  LIBDERIV_NEW_AM1=$((2*($withval)))
])
AC_SUBST(LIBDERIV_NEW_AM1)
LIBDERIV_NEW_AM2=6
AC_ARG_WITH(max-am-deriv2,[  --with-max-am-deriv2    maximum AM in ERI second derivatives],[
  LIBDERIV_NEW_AM2=$((2*($withval)))
])
AC_SUBST(LIBDERIV_NEW_AM2)

# Options for libr12
LIBR12_NEW_AM=6
LIBR12_OPT_AM=6
AC_ARG_WITH(max-am-r12,[  --with-max-am-r12       maximum AM in R12 integrals],[
  LIBR12_NEW_AM=$((2*($withval)))
  LIBR12_OPT_AM=$((2*($withval)))
])
AC_SUBST(LIBR12_NEW_AM)
AC_SUBST(LIBR12_OPT_AM)

PSI_INTEGRALS_STANDARD=1
AC_ARG_WITH(integrals-standard,[  --with-integrals-standard     specifies the standard for the basis normalization, etc.
                                The allowed values are cca (default) and psi.],[
case $withval in
  psi)
    PSI_INTEGRALS_STANDARD=0
    ;;
  cca)
    PSI_INTEGRALS_STANDARD=1
    ;;
  *)
    AC_MSG_WARN([did not recognize value $withval. will assume the default]) 
    ;;
esac
])
AC_DEFINE_UNQUOTED(PSI_INTEGRALS_STANDARD,$PSI_INTEGRALS_STANDARD)

########
## Check for programs
########

case $target_os in
  aix*)
    CC_LIST="gcc cc_r cc"
    CXX_LIST="g++ xlC_r c++"
    F77_LIST="ifort gfortran g77 xlf_r f77"
    ;;
  *)
    CC_LIST="gcc icc cc"
    CXX_LIST="g++ c++ icpc cxx"
    F77_LIST="gfortran g77 ifort f77"
    ;;
esac

AC_ARG_WITH(aix64,
[  --with-aix64            Compiling on a 64-bit AIX machine.],
CXXFLAGS="$CXXFLAGS -maix64"
F77FLAGS="$F77FLAGS -maix64"
CFLAGS="$CFLAGS -maix64"
LDFLAGS="$LDFLAGS -maix64"
AIX64=yes
AC_MSG_RESULT([Compiling on 64-bit AIX])
)

AC_ARG_WITH(cc,[  --with-cc               which C compiler to use],[
CC=$withval
])
AC_PROG_CC([$CC_LIST])

AC_ARG_WITH(cxx,[  --with-cxx              which C++ compiler to use],[
CXX=$withval
])
AC_PROG_CXX([$CXX_LIST])

AC_ARG_WITH(fc,[  --with-fc               which Fortran compiler to use],[
case $withval in
  no)
    F77="no"
    ;;
  *)
    F77=$withval
    ;;
esac
])

if test -z "$F77"; then
  AC_PROG_F77
elif test "$F77" = "no"; then
  F77=""
fi

dnl ----------- See if Fortran works --------------

if test -n "$F77" -a "$F77" != no ; then
  AC_MSG_CHECKING([if fortran compiler works])
  /bin/rm -f ffunc.f
  echo "      program main" > ffunc.f
  echo "      end" >> ffunc.f
  if $F77 -o ffunc ffunc.f 1>&5 2>&5; then
    AC_MSG_RESULT([yes])
    F77_WORKS=yes
  else
    AC_MSG_RESULT([no])
    F77_WORKS=no
  fi
  /bin/rm -f ffunc ffunc.f
fi

dnl ----------- Fortran symbol names --------------

AC_ARG_WITH(f77symbol,[  --with-f77symbol        FORTRAN77 symbol convention
                          Allowed values are:
                            lc  : lower-case
                            lcu : lower-case with underscore (default)
                            uc  : upper-case                            
                            ucu : upper-case with underscore],[
case $withval in            
  lc)
    FC_SYMBOL=1
    ;;
  lcu)
    FC_SYMBOL=2
    ;;
  uc)
    FC_SYMBOL=3
    ;;
  ucu)
    FC_SYMBOL=4
    ;;
  *)
    AC_MSG_WARN([did not recognize value $withval. will assume lower-case with underscore])    
    FC_SYMBOL=2
    ;;
esac
])

if test -z "$FC_SYMBOL"; then
 AC_MSG_CHECKING([fortran symbols])
 if test -n "$F77" -a "$F77" != no -a X$F77_WORKS != Xno; then
  /bin/rm -f ffunc.f flink.cc
  echo "      subroutine ffunc()" > ffunc.f
  echo "      return" >> ffunc.f
  echo "      end" >> ffunc.f
  $F77 -c ffunc.f 1>/dev/null 2>/dev/null
  echo "extern \"C\" void FF(); int main(int argc, char** argv){ FF(); return 0; }" > flink.cc
  if $CXX -o flink -DFF=ffunc flink.cc ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT([lowercase])
    FC_SYMBOL=1
  elif $CXX -o flink -DFF=ffunc_ flink.cc ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT([lowercase with underscore])
    FC_SYMBOL=2
  elif $CXX -o flink -DFF=FFUNC flink.cc ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT([uppercase])
    FC_SYMBOL=3
  elif $CXX -o flink -DFF=FFUNC_ flink.cc ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT([uppercase with underscore])
    FC_SYMBOL=4
  else
    AC_MSG_RESULT([giving up])
    AC_MSG_ERROR([could not determine F77 symbol names])
  fi
  /bin/rm -f ffunc.f ffunc.o flink flink.cc flink.o ffunc
 else
  FC_SYMBOL=2
  AC_MSG_RESULT([will assume lowercase with underscore])
 fi
fi
AC_SUBST(FC_SYMBOL)

AC_PROG_RANLIB
AC_PROG_CPP
AC_PROG_CXXCPP
AC_PROG_LEX([yywrap])
AC_PROG_YACC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_PROGS(AR,ar,ar) 
AC_ARG_WITH(ar,
[  --with-ar               Names the archive creator.],
AR=$withval
)
AC_SUBST(AR)

ARFLAGS=r
AC_ARG_WITH(ar-flags,
[  --with-ar-flags         Flags for the the archive creator.],
ARFLAGS=$withval
)
case $os in
  aix*)
    ARFLAGS="-X32_64 $ARFLAGS"
    ;;
esac
AC_SUBST(ARFLAGS)

AC_ARG_WITH(ld,
[  --with-ld               Names the object linker.],
LD=$withval
)
AC_SUBST(LD)

EXTRAINCDIRS=
CPPFLAGS=
AC_ARG_WITH(incdirs,
[  --with-incdirs          Specifies include directories (-Idir1 -Idir2).],
EXTRAINCDIRS=$withval
CPPFLAGS=$withval
AC_MSG_RESULT([Using extra include directories: $withval])
)

AC_ARG_WITH(libs,
[  --with-libs             Specifies libraries (-llib1 -llib2).],
LDLIBS=$withval
AC_MSG_RESULT([Using extra libraries: $withval])
)
AC_SUBST(LDLIBS)

LDFLAGS=
LIBDIRS=
AC_ARG_WITH(libdirs,
[  --with-libdirs          Specifies library directories (-Ldir1 -Ldir2).],
LIBDIRS=$withval
LDFLAGS=$withval
AC_MSG_RESULT([Using extra library directories: $withval])
)
AC_SUBST(LDFLAGS)

AC_CHECK_PROGS(PERL,perl,perl) AC_SUBST(PERL)

AC_PROG_LATEX
AC_PROG_LATEX2HTML
AC_PROG_DVIPS
AC_PROG_BIBTEX

####
# Check for doxygen and it buddies
####
AC_CHECK_PROGS(DOXYGEN,doxygen,NO)
DOT=yes
AC_ARG_WITH(dot,
[  --with-dot              Gives the path to the dot graph generator.],
DOT=$withval
)
DOT_PATH=""
if test X$DOT = Xno; then
  HAVE_DOT=NO
elif test X$DOT = Xyes; then
  AC_CHECK_PROG(HAVE_DOT,dot,YES,NO)
else
  HAVE_DOT=YES
  DOT_PATH=$DOT
fi
AC_SUBST(DOT_PATH)
AC_SUBST(HAVE_DOT)

########
## Check for libraries
########

if test ! -z "$F77"; then
  AC_F77_LIBRARY_LDFLAGS
  dnl --- This just implicitly checks that FLIBS can co-exist with the C++ linker
  AC_CHECK_LIB(c,malloc,[],[
      FLIBS=
      AC_MSG_WARN([F77 libraries don't seem to link with C++ programs])
      ],
      $FLIBS)
fi

AC_ARG_WITH(blas,[  --with-blas             which blas library to use],[
case $withval in
  yes)
    BLAS='-lblas'
    ;;
  no)
    BLAS=''
    ;;
  *)
    BLAS=$withval
    ;;
esac
])

AC_ARG_WITH(lapack,[  --with-lapack           which lapack library to use],[
case $withval in
  yes)
    LAPACK='-llapack'
    ;;
  no)
    LAPACK=''
    ;;
  *)
    LAPACK=" $withval"
    ;;
esac
])

# Set strings for particular BLAS and LAPACK functions for later testing
case $FC_SYMBOL in
  1)
    DGEMM_SYMBOL=dgemm
    DGEEV_SYMBOL=dgeev
    ;;
  2)
    DGEMM_SYMBOL=dgemm_
    DGEEV_SYMBOL=dgeev_
    ;;
  3)
    DGEMM_SYMBOL=DGEMM
    DGEEV_SYMBOL=DGEEV
    ;;
  4)
    DGEMM_SYMBOL=DGEMM_
    DGEEV_SYMBOL=DGEEV_
    ;;
esac

SAVE_LIBS=$LIBS
LIBS="$LIBS $CLIBS $FLIBS"

### If the user didn't give us a BLAS lib, see if we can find one

# If we're using the Intel compilers, we prefer MKL as well
if test -z "$BLAS"; then
  if test $CC = icc -o $CXX = icpc; then
    AC_CHECK_LIB(mkl, $DGEMM_SYMBOL, BLAS="-lmkl")
  fi
fi
# Otherwise, run through a list of possibilities
if test -z "$BLAS"; then
  AC_CHECK_LIB(essl, $DGEMM_SYMBOL, BLAS="-lessl -lpthread")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(scs, $DGEMM_SYMBOL, BLAS="-lscs")
fi

if test -z "$BLAS"; then
  SAVE_LIBS2=$LIBS
  LIBS="-lpthread $LIBS"
  AC_CHECK_LIB(goto, $DGEMM_SYMBOL, BLAS="-lgoto -lpthread")
  LIBS=$SAVE_LIBS2
fi
if test -z "$BLAS"; then
  SAVE_LIBS2=$LIBS
  LIBS="-latlas $LIBS"
  AC_CHECK_LIB(f77blas, $DGEMM_SYMBOL, BLAS="-lf77blas -latlas")
  LIBS=$SAVE_LIBS2
fi  
if test -z "$BLAS"; then
  AC_CHECK_LIB(mkl, $DGEMM_SYMBOL, BLAS="-lmkl")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(blas, $DGEMM_SYMBOL, BLAS="-lblas")
fi
if test -z "$BLAS"; then
case $target_vendor in
  apple)
    BLAS="-altivec -framework vecLib"
    AUTOAPPLE="yes"
    AC_MSG_RESULT([No BLAS specified and vendor is Apple: Assuming vecLib is available!])
    ;;
esac
fi
if test -z "$BLAS"; then
  AC_MSG_ERROR("Did not find a BLAS library.")
fi
AC_SUBST(BLAS)
AC_MSG_RESULT([BLAS library is $BLAS])
LIBS=$SAVE_LIBS

# Similar checks for LAPACK, except that sometimes BLAS includes LAPACK
# Also note that LAPACK often (always?) depends on BLAS
SAVE_LIBS=$LIBS
LIBS="$BLAS $LIBS $CLIBS $FLIBS"
if test -z "$LAPACK"; then
  if test ! -z "$AUTOAPPLE"; then
    AC_MSG_RESULT([Assuming vecLib for LAPACK])
    LAPACK="$BLAS"
    LAPACK_IN_BLAS=yes
  else
    AC_CHECK_FUNC($DGEEV_SYMBOL,LAPACK_IN_BLAS=yes)
  fi
fi 
if test -z "$LAPACK_IN_BLAS"; then
  if test -z "$LAPACK"; then
    if test $CC = icc -o $CXX = icpc; then
      # Intel's libmkl_lapack64 contains the double-precision routines
      AC_CHECK_LIB(mkl_lapack64,$DGEEV_SYMBOL,LAPACK="-lmkl_lapack64")
    fi
  fi
  if test -z "$LAPACK"; then
    AC_CHECK_LIB(lapack,$DGEEV_SYMBOL,LAPACK="-llapack")
  fi
  if test -z "$LAPACK"; then
    AC_MSG_ERROR("Did not find a LAPACK library.");
  fi
  AC_MSG_RESULT([LAPACK library is $LAPACK])
else
  AC_MSG_RESULT([LAPACK library is included in BLAS+LIBS])
fi
AC_SUBST(LAPACK)
LIBS=$SAVE_LIBS

# Check for libcompat.a availability
#case $target_vendor in
#  apple)
#    AC_CHECK_LIB([compat],[main],[CLIBS="-lcompat $CLIBS"],[AC_MSG_ERROR(Missing libcompat.a: See http://www.opensource.apple.com/ for more information)],[])ac_cv_lib_compat=ac_cv_lib_compat_main

#    echo "$CLIBS"
#  ;;
#esac

########
## Check for header files
########

# check for header files
AC_CHECK_HEADERS([errno.h fcntl.h unistd.h stdint.h sys/types.h sys/times.h sys/param.h])
# check for C++ headers
AC_CXX_PROCESS_CHECK([cmath],[#include <cmath>],,[HAVE_CMATH=1],[HAVE_CMATH=0]) 
#CXXDEF="$CXXDEF -DHAVE_CMATH=$HAVE_CMATH"
AC_DEFINE_UNQUOTED(HAVE_CMATH,$HAVE_CMATH)
AC_CXX_PROCESS_CHECK([cstdio],[#include <cstdio>],,[HAVE_CSTDIO=1],[HAVE_CSTDIO=0]) 
#CXXDEF="$CXXDEF -DHAVE_CSTDIO=$HAVE_CSTDIO"
AC_DEFINE_UNQUOTED(HAVE_CSTDIO,$HAVE_CSTDIO)
AC_CXX_PROCESS_CHECK([cstdlib],[#include <cstdlib>],,[HAVE_CSTDLIB=1],[HAVE_CSTDLIB=0]) 
#CXXDEF="$CXXDEF -DHAVE_CSTDLIB=$HAVE_CSTDLIB"
AC_DEFINE_UNQUOTED(HAVE_CSTDLIB,$HAVE_CSTDLIB)
AC_CXX_PROCESS_CHECK([cstring],[#include <cstring>],,[HAVE_CSTRING=1],[HAVE_CSTRING=0]) 
#CXXDEF="$CXXDEF -DHAVE_CSTRING=$HAVE_CSTRING"
AC_DEFINE_UNQUOTED(HAVE_CSTRING,$HAVE_CSTRING)
AC_CXX_PROCESS_CHECK([cstddef],[#include <cstddef>],,[HAVE_CSTDDEF=1],[HAVE_CSTDDEF=0]) 
#CXXDEF="$CXXDEF -DHAVE_CSTDDEF=$HAVE_CSTDDEF"
AC_DEFINE_UNQUOTED(HAVE_CSTDDEF,$HAVE_CSTDDEF)

####################
#  Check for Ruby  #
####################
HAVE_RUBY="no"
RUBYINC=""
CHECKRUBY="yes"
AC_ARG_WITH(ruby-include,[  --with-ruby-include     where to look for Ruby headers (-Idir1 -Idir2 etc.)],[
case $withval in
  yes)
    ;;
  no)
        CHECKRUBY="no"
    ;;
  *)
        RUBYINC="$withval"
    ;;
esac
])
if test X$CHECKRUBY = Xyes; then
        REF_CPPFLAGS=$CPPFLAGS CPPFLAGS="$REF_CPPFLAGS $RUBYINC"
        AC_CHECK_HEADER(ruby.h,[HAVE_RUBY="yes"],[HAVE_RUBY="no"])
        CPPFLAGS=$REF_CPPFLAGS
fi

if test X$HAVE_RUBY = Xyes; then
  HAVE_RUBY="no"
  RUBYLIB="-lruby-static"
  CHECKRUBY="yes"
  AC_ARG_WITH(ruby-lib,[  --with-ruby-lib         which Ruby library to use (-Ldir1 -lruby_red)],[
  case $withval in
  yes)
    ;;
  no)
	CHECKRUBY="no"
    ;;
  *)
        RUBYLIB="$withval"
    ;;
  esac
  ])
  if test $CHECKRUBY = yes; then
  AC_MSG_CHECKING([for ruby_init])
  REF_LIBS=$LIBS  LIBS="$RUBYLIB $LIBS"
  REF_CPPFLAGS=$CPPFLAGS CPPFLAGS="$REF_CPPFLAGS $RUBYINC"
  AC_LINK_IFELSE([
    AC_LANG_PROGRAM(
      [[#include <ruby.h>
      ]],
      [[
        void (*funcptr)() = ruby_init;
        funcptr();
      ]]
      ),
    HAVE_RUBY=yes
    AC_MSG_RESULT(yes),
    HAVE_RUBY=no
    AC_MSG_RESULT(no)
  ])
  LIBS=$REF_LIBS
  CPPFLAGS=$REF_CPPFLAGS
  fi
fi

#########################################################################
# Special linking flags for psirb, because it loads shared object files #
#  This has to do with dlopen being "lazy" by default on Apple, but not #
#  on Linux                                                             #
#########################################################################
case $target_vendor in
	apple)
		RUBYLDFLAGS=""
		;;
	*)
		RUBYLDFLAGS="-rdynamic -Wl,-export-dynamic"
		;;
esac
AC_SUBST(RUBYLIB)
AC_SUBST(RUBYINC)
AC_SUBST(HAVE_RUBY)
AC_DEFINE_UNQUOTED(HAVE_RUBY,$HAVE_RUBY)
AC_SUBST(RUBYLDFLAGS)
########
## Check for compiler characteristics
########

test -z "$LIBSUF" && LIBSUF=a AC_SUBST(LIBSUF)
test -z "$OBJSUF" && OBJSUF=o AC_SUBST(OBJSUF)

dnl -- check how dependency information is built --
# The GNU compilers work with:
CCDEPENDSUF=none
CXXDEPENDSUF=none
CCDEPENDFLAGS=-M
CXXDEPENDFLAGS=-M
/bin/rm -f depcheck.u depcheck.c depcheck.cc depcheck.o
# Check for an IBM visual age C compiler
echo "#include <stdio.h>" > depcheck.c
$CC $CPPFLAGS $CFLAGS -M depcheck.c > /dev/null 2>&1
if test -f depcheck.u; then
  CCDEPENDSUF=u
fi
/bin/rm -f depcheck.u depcheck.c depcheck.cc depcheck.o
# Check for an IBM visual age C++ compiler
echo "#include <iostream>" > depcheck.cc
$CXX $CPPFLAGS $CXXFLAGS -M -E depcheck.cc > /dev/null 2>&1
if test -f depcheck.u; then
  CXXDEPENDSUF=u
  CXXDEPENDFLAGS="-M -E"
fi
/bin/rm -f depcheck.u depcheck.c depcheck.cc depcheck.o

AC_SUBST(CCDEPENDSUF)
AC_SUBST(CXXDEPENDSUF)
AC_SUBST(CCDEPENDFLAGS)
AC_SUBST(CXXDEPENDFLAGS)

AC_CHECK_DECLS([putenv],[extradef="-DHAVE_DECL_PUTENV=1"])
#CDEF="$CDEF $extradef"
#CXXDEF="$CXXDEF $extradef"
AC_DEFINE_UNQUOTED(HAVE_DECL_PUTENV,1)

AC_CHECK_FUNC(putenv,[extradef="-DHAVE_PUTENV=1"])
#CDEF="$CDEF $extradef"
#CXXDEF="$CXXDEF $extradef"
AC_DEFINE_UNQUOTED(HAVE_PUTENV,1)

AC_CHECK_DECLS([setenv],[extradef="-DHAVE_DECL_SETENV=1"])
#CDEF="$CDEF $extradef"
#CXXDEF="$CXXDEF $extradef"
AC_DEFINE_UNQUOTED(HAVE_DECL_SETENV,1)

AC_CHECK_FUNC(setenv,[extradef="-DHAVE_SETENV=1"])
#CDEF="$CDEF $extradef"
#CXXDEF="$CXXDEF $extradef"
AC_DEFINE_UNQUOTED(HAVE_SETENV,1)

#AC_CHECK_LIB(m,isinf,[extradef="-DHAVE_FUNC_ISINF=1"],[extradef=""])
AC_CHECK_LIB(m,isinf,[AC_DEFINE_UNQUOTED(HAVE_FUNC_ISINF,1)],[
  AC_CHECK_LIB(m,finite,[AC_DEFINE_UNQUOTED(HAVE_FUNC_FINITE,1)],[
    AC_MSG_ERROR([Did not find isinf() or finite() functions.])
  ])
])
#if test X$extradef = X; then
#  AC_CHECK_LIB(m,finite,[AC_DEFINE_UNQUOTED(HAVE_FUNC_FINITE,1)],[extradef=""])
#fi
#if test X$extradef = X; then
#  AC_MSG_ERROR([Did not find isinf() or finite() functions.])
#else
  #CDEF="$CDEF $extradef"
  #CXXDEF="$CXXDEF $extradef"
#fi


dnl ------------------- Determine Compiler Flags --------------------

# User's input overrides all defaults
# Defaults depend on the compiler and the architecture
AC_MSG_RESULT([-------------------------------------------------------------------------------])

# default is to optimize with -O2, no debugging symbols will be generated
PSICPPFLAGS='-I. -I$(objincdir) -I$(top_srcdir)/include -I$(top_srcdir)/src/lib -I$(top_objdir)/src/lib' AC_SUBST(PSICPPFLAGS)
CDBG='' AC_SUBST(CDBG)
COPT='-O2' AC_SUBST(COPT)
AC_SUBST(CDEF)
COTH='' AC_SUBST(COTH)
CINC='$(CPPFLAGS) $(PSICPPFLAGS)' AC_SUBST(CINC)
CXXDBG='' AC_SUBST(CXXDBG)
CXXOPT='-O2' AC_SUBST(CXXOPT)
AC_SUBST(CXXDEF)
CXXOTH='' AC_SUBST(CXXOTH)
CXXINC='$(CPPFLAGS) $(PSICPPFLAGS)' AC_SUBST(CXXINC)

# Check user's input
# The default is to optimize.  --with-opt=no or --without-opt will turn
# off optimization.
OPTFLAGS_GIVEN=no
AC_ARG_WITH(opt,[  --with-opt              Specifies the C++/C optimization compiler flags.
                          --with-opt=no means no optimization],[
case $withval in
  no)
    COPT=''
    CXXOPT=''
    OPTFLAGS_GIVEN=yes
    ;;
  *)
    COPT="$withval"
    CXXOPT="$withval"
    OPTFLAGS_GIVEN=yes
    ;;
esac
])
if test X$OPTFLAGS_GIVEN = Xyes; then
  AC_MSG_RESULT([C++ optimization flags           : $CXXOPT])
fi

# if not given by the user, use architecture-dependent defaults
if test X$OPTFLAGS_GIVEN = Xno; then
 case $target_cpu in
  rs6000)
    if test "$CC" = xlc_r -o "$CC" = cc_r ; then
      COPT="-O3 -qhot -qstrict -qarch=auto -qtune=auto -qmaxmem=-1"
    fi
    if test "$CXX" = xlC_r ; then
      CXXOPT="-O3 -qhot -qstrict -qarch=auto -qtune=auto -qmaxmem=-1"
    fi
    ;;
  powerpc)
    if test "$CC" = xlc_r -o "$CC" = cc_r ; then
      COPT="-O3 -qhot -qstrict -qarch=auto -qtune=auto -qmaxmem=-1"
    fi
    if test "$CXX" = xlC_r ; then
      CXXOPT="-O3 -qhot -qstrict -qarch=auto -qtune=auto -qmaxmem=-1"
      LDLIBS="$LDLIBS $FLIBS"
    fi
    if test "$CXX" = g++ ; then
      CXXOPT="-O3"
      LDLIBS="$LDLIBS $FLIBS"
    fi
    ;;
  x86_64)
    if test "$CC" = gcc -o "$CC" = cc ; then
      COPT="-O3 -m64"
    fi
    if test "$CXX" = g++ -o "$CXX" = c++ -o "$CXX" = gcc ; then
      CXXOPT="-O3 -m64"
    fi
    ;;
 esac
 AC_MSG_RESULT([C++ optimization flags           : $CXXOPT])
fi

case $target_os in
  aix*)
# _ALL_SOURCE is required because some include files on AIX don't seem to work with pure POSIX
# Removed -D_LARGE_FILE 
    EXTRADEFS="-D_ALL_SOURCE"
    CDEF="$CDEF $EXTRADEFS"
    CXXDEF="$CXXDEF $EXTRADEFS"
    ;;
  linux*)
    EXTRADEFS="-DLinux -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE"
    CDEF="$CDEF $EXTRADEFS"
    CXXDEF="$CXXDEF $EXTRADEFS"
    ;;
esac
AC_MSG_RESULT([extra preprocessor defines       : $CXXDEF])

# Obtain the debug flags.
# The default is to not debug.  --with-debug will turn on debugging.
AC_ARG_WITH(debug,[  --with-debug            will turn on debugging],[
case $withval in
  yes)
    CDBG='-g'
    CXXDBG='-g'
    ;;
  no)
    CDBG=''
    CXXDBG=''
    ;;
  *)
    CDBG="$withval"
    CXXDBG="$withval"
    ;;
esac
])
AC_MSG_RESULT([debugging C++ compiler flags     : $CXXDBG])

# Adjust the flags dealing with strictness of the compilers.
# The default is to disable strict language checks.  --with-opt=no or 
# --without-opt will turn off optimization.
AC_ARG_WITH(strict,[  --with-strict           --with-strict=yes to turn on strict compiler warnings],[
case $withval in
  yes)
    # "Strict" flags only known for GNU compilers
    if test "$GCC" = yes ; then
      COTH="$COTH -Wall -ansi"
    fi
    if test "$GXX" = yes ; then
      CXXOTH="$CXXOTH -Wall -ansi"
    fi
    ;;
  no)
    ;;
esac
])
AC_MSG_RESULT([misc C++ compiler flags          : $CXXOTH])

CFLAGS='$(CDBG) $(COPT) $(CDEF) $(COTH) $(CINC)' AC_SUBST(CFLAGS)
CXXFLAGS='$(CXXDBG) $(CXXOPT) $(CXXDEF) $(CXXOTH) $(CXXINC)' AC_SUBST(CXXFLAGS)

if test X$AIX64 = Xyes; then
  if test X$CXX = Xg++; then
    COTH="$COTH -maix64" 
    CXXOTH="$CXXOTH -maix64" 
    F77FLAGS="$F77FLAGS -maix64" 
    LDFLAGS="$LDFLAGS -maix64" 
    LDLIBS="$LAPACK $BLAS $LDLIBS"
  fi
fi

AC_MSG_RESULT([-------------------------------------------------------------------------------])


########
## Output
########

AC_CONFIG_FILES([include/Makefile
          include/chkpt_params.h
          Makefile lib/Makefile lib/ruby/Makefile src/Makefile
          src/lib/Makefile src/lib/MakeRules src/lib/MakeVars
          src/util/Makefile src/util/MakeRules src/util/MakeVars
          src/bin/Makefile src/bin/MakeRules src/bin/MakeVars
          doc/MakeVars doc/MakeRules doc/Makefile
          doc/macros.tex doc/progman/svn.eps
          doc/installation/Makefile
          doc/progman/Makefile doc/userman/Makefile
          doc/userman/prsty.bst
          doc/doxygen/Makefile
          doc/doxygen/doxygen.cfg

          src/lib/libdpd/Makefile
          src/lib/libiwl/Makefile
          src/lib/libpsio/Makefile
          src/lib/libqt/Makefile
          src/lib/libbasis/Makefile
          src/lib/libciomr/Makefile
          src/lib/libipv1/Makefile
          src/lib/libchkpt/Makefile
          src/lib/libderiv/Makefile
          src/lib/libderiv/MakeVars
          src/lib/libderiv/libderiv_config.h
          src/lib/libint/Makefile
          src/lib/libint/MakeVars
          src/lib/libint/libint_config.h
          src/lib/libr12/Makefile
          src/lib/libr12/MakeVars
          src/lib/libr12/libr12_config.h
          src/lib/libutil/Makefile
          src/lib/liboptions/Makefile
          src/lib/libmoinfo/Makefile
          src/lib/libmints/Makefile
          src/util/tocprint/Makefile
          src/util/psi2molden/Makefile

          src/bin/cphf/Makefile
          src/bin/psi3/Makefile
          src/bin/ccdensity/Makefile
          src/bin/ccenergy/Makefile
          src/bin/cchbar/Makefile
          src/bin/cclambda/Makefile
          src/bin/ccsort/Makefile
          src/bin/cctriples/Makefile
          src/bin/cceom/Makefile
          src/bin/ccresponse/Makefile
          src/bin/cints/Makefile
          src/bin/cis/Makefile
          src/bin/clag/Makefile
          src/bin/cscf/Makefile
          src/bin/detci/Makefile
          src/bin/detcas/Makefile
          src/bin/detcasman/Makefile
          src/bin/extrema/Makefile
          src/bin/geom/Makefile
          src/bin/input/Makefile
          src/bin/intder/Makefile
          src/bin/mcscf/Makefile
          src/bin/mp2/Makefile
          src/bin/mp2r12/Makefile
          src/bin/oeprop/Makefile
          src/bin/optking/Makefile
          src/bin/psiclean/Makefile
          src/bin/transqt/Makefile
          src/bin/transqt2/Makefile
          src/bin/mvo/Makefile
          src/bin/dboc/Makefile
          src/bin/cusp/Makefile
          src/bin/mocube/Makefile
          src/bin/localize/Makefile
          src/bin/stable/Makefile
          src/bin/response/Makefile
          src/bin/psirb/Makefile
          src/bin/psimrcc/Makefile
          src/bin/nonbonded/Makefile
          src/bin/thermo/Makefile
          src/bin/df-mp2/Makefile
          src/bin/scf/Makefile

          tests/Makefile
          tests/MakeVars
          tests/MakeRules
          tests/runtest.pl
          tests/scf-opt/Makefile
          tests/scf-opt-fragments/Makefile
          tests/mp2-opt-fragments/Makefile
          tests/scf-opt2/Makefile
          tests/scf-opt3/Makefile
          tests/scf-opt4/Makefile
          tests/scf-opt5/Makefile
          tests/scf-opt6/Makefile
          tests/scf-opt7/Makefile
          tests/scf-opt9/Makefile
          tests/scf-opt10/Makefile
          tests/scf-opt11/Makefile
          tests/scf-opt12/Makefile
          tests/scf-opt13/Makefile
          tests/scf-opt14/Makefile
          tests/scf-opt-numer/Makefile
          tests/scf-opt2-numer/Makefile
          tests/scf-freq/Makefile
          tests/scf-freq-first-3/Makefile
          tests/scf-freq-first-3a/Makefile
          tests/scf-freq-first-5/Makefile
          tests/scf-freq-first-5a/Makefile
          tests/scf-freq-none-3/Makefile
          tests/scf-freq-none-5/Makefile
          tests/scf-symm-fc-numer/Makefile
          tests/scf-fc-numer/Makefile
          tests/scf-polar/Makefile
          tests/rhf-oeprop/Makefile
          tests/rhf-stab/Makefile
          tests/rhf-stab/rhf-stab.pl
          tests/uhf-stab/Makefile
          tests/uhf-stab/uhf-stab.pl
          tests/rohf-stab/Makefile
          tests/rohf-stab/rohf-stab.pl
          tests/cis-sp/Makefile
          tests/cisd-h2o+-0/Makefile
          tests/cisd-h2o+-1/Makefile
          tests/cisd-h2o+-2/Makefile
          tests/cisd-sp/Makefile
          tests/cisd-sp-2/Makefile
          tests/cisd-opt-numer/Makefile
          tests/fci-h2o/Makefile
          tests/fci-h2o-2/Makefile
          tests/fci-h2o-clpse/Makefile
          tests/fci-h2o-fzcv/Makefile
          tests/fci-dipole/Makefile
	  tests/fci-trans/Makefile
          tests/casscf-sp/Makefile
          tests/casscf-sa-sp/Makefile
          tests/casscf-fzc-sp/Makefile
          tests/rasscf-sp/Makefile
          tests/casscf-opt-numer/Makefile
          tests/casscf-opt/Makefile
          tests/zaptn-nh2/Makefile
          tests/mp2-sp/Makefile
          tests/mp2-direct-sp/Makefile
          tests/mp2-oeprop/Makefile
          tests/mp2-scs/Makefile
          tests/cc1/Makefile
          tests/cc2/Makefile
          tests/cc3/Makefile
          tests/cc4/Makefile
          tests/cc5/Makefile
          tests/cc6/Makefile
          tests/cc7/Makefile
          tests/cc8/Makefile
          tests/cc8a/Makefile
          tests/cc8b/Makefile
          tests/cc8c/Makefile
          tests/cc9/Makefile
          tests/cc9a/Makefile
          tests/cc10/Makefile
          tests/cc11/Makefile
          tests/cc12/Makefile
          tests/cc13/Makefile
          tests/cc13a/Makefile
          tests/cc14/Makefile
          tests/cc15/Makefile
          tests/cc16/Makefile
          tests/cc17/Makefile
          tests/cc18/Makefile
          tests/cc19/Makefile
          tests/cc21/Makefile
          tests/cc22/Makefile
          tests/cc23/Makefile
          tests/cc24/Makefile
          tests/cc25/Makefile
          tests/cc26/Makefile
          tests/cc27/Makefile
          tests/cc28/Makefile
          tests/cc29/Makefile
          tests/cc30/Makefile
          tests/cc31/Makefile
          tests/cc32/Makefile
          tests/cc33/Makefile
          tests/cc34/Makefile
          tests/cc35/Makefile
          tests/cc36/Makefile
          tests/cc37/Makefile
          tests/cc38/Makefile
          tests/cc39/Makefile
          tests/cc40/Makefile
          tests/cc41/Makefile
          tests/cc42/Makefile
          tests/cc43/Makefile
          tests/cc44/Makefile
          tests/cc45/Makefile
          tests/cc46/Makefile
          tests/cc47/Makefile
          tests/cc48/Makefile
          tests/cc49/Makefile
          tests/cc50/Makefile
          tests/cc51/Makefile
          tests/cc52/Makefile
          tests/ooccd/Makefile
          tests/pccsd1/Makefile
          tests/psimrcc-sp1/Makefile
          tests/psimrcc-idmrpt2/Makefile
          tests/psimrcc-mp2-ccsd-1/Makefile
          tests/psimrcc-mp2-ccsd-2/Makefile
          tests/psimrcc-mp2-ccsd-3/Makefile
          tests/psimrcc-mp2-ccsd-4/Makefile
          tests/psimrcc-ccsd_t-1/Makefile
          tests/psimrcc-ccsd_t-2/Makefile
          tests/psimrcc-ccsd_t-3/Makefile
          tests/psimrcc-ccsd_t-4/Makefile
          tests/dboc-rhf1/Makefile
          tests/dboc-rohf1/Makefile
          tests/dboc-uhf1/Makefile
          tests/dboc-rcisd1/Makefile
          tests/dboc-rocisd1/Makefile
          tests/extrema-zmat/Makefile
          tests/extrema-deloc/Makefile
          tests/mcscf-rhf1/Makefile
          tests/mcscf-rohf1/Makefile
          tests/mcscf-twocon1/Makefile
          tests/mp2r12-sp1/Makefile
          tests/mp2-opt-numer/Makefile
          tests/mp2-opt/Makefile
          tests/rhf-lindep1/Makefile
          tests/rohf-lindep1/Makefile
          tests/uhf-lindep1/Makefile
          tests/psi_start/Makefile
          tests/psi_start/psi_start.pl
          tests/rb-scf-opt-numer/Makefile
          tests/rb-ccsd_t-sp/Makefile
          tests/rb-fci-h2o/Makefile
          tests/rb-ccsd-opt-grad/Makefile
          tests/scf-mvd-opt/Makefile
          tests/scf-mvd-opt-numer/Makefile
          tests/scf-mvd-opt-puream/Makefile
          tests/scf-mvd-opt-numer-puream/Makefile
          tests/scf+d-opt1/Makefile
          tests/scf+d-freq/Makefile
          tests/rhf-efield/Makefile
          tests/casscf-sp-reorder/Makefile

          src/bin/giao/Makefile

          src/samples/Makefile
          src/samples/MakeVars
          src/samples/MakeRules
          src/samples/mints/Makefile
])

########################
#  Check for loopmrcc  #
########################

COMPILE_LOOPMRCC="no"
AC_ARG_WITH(loopmrcc,[  --with-loopmrcc         compile loopmrcc],[
case $withval in
  yes)
        COMPILE_LOOPMRCC="yes"
    ;;
  no)
    ;;
  *)
    ;;
esac
])

if test X$COMPILE_LOOPMRCC = Xyes; then
        AC_CONFIG_FILES([src/bin/loopmrcc_old/Makefile])
        AC_CONFIG_FILES([src/bin/loopmrcc/Makefile])
fi

AC_SUBST(COMPILE_LOOPMRCC)

COMPILE_WICKED="no"
AC_ARG_WITH(wicked,[  --with-wicked           compile wick&d],[
case $withval in
  yes)
        COMPILE_WICKED="yes"
    ;;
  no)
    ;;
  *)
    ;;
esac
])

if test X$COMPILE_WICKED = Xyes; then
        AC_CONFIG_FILES([src/bin/wick&d/Makefile])
fi

AC_SUBST(COMPILE_WICKED)

AC_OUTPUT

if test $HAVE_RUBY = no; then
	if test -z $RUBYLIB; then
		AC_MSG_NOTICE([Unable to locate libruby-static.a will not link the psirb module.])
	fi
	if test -z $RUBYINC; then
		AC_MSG_NOTICE([Unable to locate ruby.h will not link the psirb module.])
	fi
fi
