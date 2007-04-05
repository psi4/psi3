# Process this file with autoconf to produce a configure script.

define([psi_version],[3.3.0])
define([psi_bugreport],[psicode@users.sourceforge.net])
AC_PREREQ(2.57)
AC_INIT([psi],[psi_version],[psi_bugreport])

PSI_VERSION=psi_version
PSI_BUGREPORT=psi_bugreport
AC_SUBST(PSI_VERSION)
AC_SUBST(PSI_BUGREPORT)

AC_CANONICAL_TARGET
AC_PREFIX_DEFAULT(/usr/local/psi)

objincdir='$(top_objdir)/include'
AC_SUBST(objincdir)

docdir='$(prefix)/doc'
mandir='$(docdir)/man'
txtdir='$(docdir)/txt'
htmldir='$(docdir)/html'
dvidir='$(docdir)/dvi'

AC_SUBST(docdir)
AC_SUBST(mandir)
AC_SUBST(htmldir)
AC_SUBST(dvidir)
AC_SUBST(txtdir)

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
LIBDERIV_NEW_AM2=4
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

########
## Check for programs
########

case $target_os in
  aix*)
    CC_LIST="cc_r gcc cc"
    CXX_LIST="xlC_r g++ c++"
    F77_LIST="xlf_r gfortran g77 ifort f77 fc f2c"
    ;;
  irix*)
    CC_LIST="cc gcc"
    CXX_LIST="CC g++ c++"
    F77_LIST="f77 gfortran g77 ifort fc f2c"
    ;;
  *)
    CC_LIST="gcc icc cc"
    CXX_LIST="g++ c++ icpc cxx"
    F77_LIST="gfortran g77 ifort f77 fc f2c"
    ;;
esac


AC_ARG_WITH(cc,[  --with-cc               which C compiler to use],[
case $withval in
  yes)
    CC='gcc'
    ;;
  no)
    CC=''
    ;;
  *)
    CC=$withval
    ;;
esac
])

if test -z "$CC"; then
  AC_PROG_CC([$CC_LIST])
fi
AC_SUBST(CC)

AC_ARG_WITH(cxx,[  --with-cxx              which C++ compiler to use],[
case $withval in           
  yes)
    CXX='g++'
    ;;
  no)
    CXX=''
    ;;
  *)
    CXX=$withval
    ;;
esac
])

if test -z "$CXX"; then
  AC_PROG_CXX([$CXX_LIST])
fi
AC_SUBST(CXX)

AC_ARG_WITH(fc,[  --with-fc               which Fortran compiler to use],[
case $withval in
  yes)
    F77='g77'
    ;;
  no)
    F77='no'
    ;;
  *)
    F77=$withval
    ;;
esac
])

if test "$F77" != no; then
  if test -z "$F77"; then
    AC_PROG_F77([$F77_LIST])
  fi
else
  F77=
fi
AC_SUBST(FC,[$F77])

AC_PROG_RANLIB
AC_CHECK_PROGS(AR,ar,ar) AC_SUBST(AR)
AC_PROG_CPP
AC_PROG_CXXCPP
AC_PROG_YACC
AC_PROG_INSTALL
AC_PROG_LN_S

AC_CHECK_PROGS(PERL,perl,perl) AC_SUBST(PERL)


########
## Check for libraries
########

if test ! -z "$F77"; then
  AC_F77_LIBRARY_LDFLAGS
fi
# This just implicitly checks that FLIBS can co-exist with the C linker
AC_CHECK_LIB(c,exit,[],[
FLIBS=
AC_MSG_WARN([F77 libraries don't seem to link with C programs])
],$FLIBS)

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

########
## Check for header files
########

########
## Check for types
########

########
## Check for compiler characteristics
########

test -z "$F77SUF" && F77SUF=f AC_SUBST(F77SUF)
test -z "$LIBSUF" && LIBSUF=a AC_SUBST(LIBSUF)
test -z "$OBJSUF" && OBJSUF=o AC_SUBST(OBJSUF)

AC_CC_PROCESS_CHECK(C dependency processing,[
#include "nosuchincludefile.h"
],-M,CCDEPEND='$(CC) -M'; DODEPEND=yes,CCDEPEND=":"; DODEPEND=no)
AC_SUBST(CCDEPEND)
AC_SUBST(DODEPEND)
AC_MSG_RESULT(CCDEPEND is $CCDEPEND)
AC_MSG_RESULT(DODEPEND is $DODEPEND)

AC_CHECK_LIB(m,exit,CLIBS="$CLIBS -lm")
AC_CHECK_LIB(bsd,exit,CLIBS="$CLIBS -lbsd")
AC_MSG_RESULT(c libs: $CLIBS)
AC_SUBST(CLIBS)

AC_CHECK_LIB(m,exit,CXXLIBS="$CXXLIBS -lm")
AC_MSG_RESULT(c++ libs: $CXXLIBS)
AC_SUBST(CXXLIBS)


## Fortran-C linking convention tests

AC_ARG_WITH(f77symbol,[  --with-f77symbol        FORTRAN77 symbol convention (if no Fortran compiler). Allowed values are:
                            lc  (lower-case)
                            lcu (lower-case with underscore; default, if
not given)
                            uc  (upper-case)
                            ucu (upper-case with underscore)],[
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
  AC_MSG_CHECKING(Fortran symbols)
if test -n "$F77" ; then
  /bin/rm -f ffunc.f flink.c
  echo "      subroutine ffunc()" > ffunc.f
  echo "      return" >> ffunc.f
  echo "      end" >> ffunc.f
  $F77 -c ffunc.f 1>/dev/null 2>/dev/null
  echo "main(){ FF(); return 0; }" > flink.c
  if $CC -o flink -DFF=ffunc flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(same as C)
    FC_SYMBOL=1
  elif $CC -o flink -DFF=ffunc_ flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(lowercase with underscore)
    FC_SYMBOL=2
  elif $CC -o flink -DFF=FFUNC flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(uppercase)
    FC_SYMBOL=3
  elif $CC -o flink -DFF=FFUNC_ flink.c ffunc.o $LDFLAGS $LIBS 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT(uppercase with underscore)
    FC_SYMBOL=4
  else
    AC_MSG_RESULT(giving up)
    AC_MSG_ERROR(could not determine fortran symbol names)
  fi
  /bin/rm -f ffunc.f ffunc.o flink flink.c flink.o ffunc
else
  FC_SYMBOL=2
  AC_MSG_RESULT(No Fortran compiler provided; guessing lowercase with underscore)
fi
fi
AC_SUBST(FC_SYMBOL)

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

########
## Check for library functions
########

SAVE_LIBS=$LIBS
LIBS="$LIBS $CLIBS $FLIBS"
echo "LIBS = $LIBS"

# If the user didn't give us a BLAS lib, see if we can find one
if test -z "$BLAS"; then
  AC_CHECK_LIB(essl, $DGEMM_SYMBOL, BLAS="-lessl")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(scs, $DGEMM_SYMBOL, BLAS="-lscs")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(goto, $DGEMM_SYMBOL, BLAS="-lgoto")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(f77blas, $DGEMM_SYMBOL, BLAS="-lf77blas -latlas")
fi
if test -z "$BLAS"; then
  AC_CHECK_LIB(blas, $DGEMM_SYMBOL, BLAS="-lblas")
fi
if test -z "$BLAS"; then
  AC_MSG_ERROR("Did not find a BLAS library.")
fi
AC_SUBST(BLAS)
AC_MSG_RESULT([BLAS library is $BLAS])

LIBS=$SAVE_LIBS

# Similar checks for LAPACK, except that sometimes BLAS includes LAPACK
if test -z "$LAPACK"; then
  CURRENT_LIBS=$LIBS
  LIBS="$BLAS $LIBS"
  AC_CHECK_FUNC($DGEEV_SYMBOL,LAPACK_IN_BLAS=yes)
  LIBS=$CURRENT_LIBS
fi
if test -z "$LAPACK_IN_BLAS"; then
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

########
## Output
########

AC_CONFIG_FILES([include/sizeof.h include/Makefile
          include/chkpt_params.h
          Makefile lib/Makefile src/Makefile
          src/lib/Makefile src/lib/MakeRules src/lib/MakeVars
          src/util/Makefile src/util/MakeRules src/util/MakeVars
          src/bin/Makefile src/bin/MakeRules src/bin/MakeVars
          doc/MakeVars doc/MakeRules doc/Makefile
          doc/macros.tex
          doc/installation/Makefile
          doc/progman/Makefile doc/userman/Makefile

          src/lib/libderiv/Makefile
          src/lib/libderiv/MakeVars
          src/lib/libderiv/libderiv_config.h
          src/lib/libint/Makefile
          src/lib/libint/MakeVars
          src/lib/libint/libint_config.h
          src/lib/libr12/Makefile
          src/lib/libr12/MakeVars
          src/lib/libr12/libr12_config.h
])
          
#          src/util/tocprint/Makefile
#          src/util/psi2molden/Makefile
#
#          src/lib/libciomr/Makefile
#          src/lib/libciomr/iomrparam.h
#          src/lib/libipv1/Makefile
#          src/lib/libchkpt/Makefile
#
#          src/lib/libdpd/Makefile
#          src/lib/libiwl/Makefile
#          src/lib/libpsio/Makefile
#          src/lib/libqt/Makefile
#          src/lib/libbasis/Makefile
#
#          src/bin/cphf/Makefile
#          src/bin/psi3/Makefile
#          src/bin/ccdensity/Makefile
#          src/bin/ccenergy/Makefile
#          src/bin/cchbar/Makefile
#          src/bin/cclambda/Makefile
#          src/bin/ccsort/Makefile
#          src/bin/cctriples/Makefile
#          src/bin/cceom/Makefile
#          src/bin/ccresponse/Makefile
#          src/bin/cints/Makefile
#          src/bin/cis/Makefile
#          src/bin/clag/Makefile
#          src/bin/cscf/Makefile
#          src/bin/detci/Makefile
#          src/bin/detcas/Makefile
#          src/bin/detcasman/Makefile
#          src/bin/extrema/Makefile
#          src/bin/geom/Makefile
#          src/bin/input/Makefile
#          src/bin/intder/Makefile
#          src/bin/mp2/Makefile
#          src/bin/mp2r12/Makefile
#          src/bin/oeprop/Makefile
#          src/bin/optking/Makefile
#          src/bin/psiclean/Makefile
#          src/bin/transqt/Makefile
#          src/bin/transqt2/Makefile
#          src/bin/mvo/Makefile
#          src/bin/dboc/Makefile
#          src/bin/cusp/Makefile
#          src/bin/mocube/Makefile
#          src/bin/localize/Makefile
#          src/bin/stable/Makefile
#          src/bin/response/Makefile
#          tests/Makefile
#          tests/MakeVars
#          tests/MakeRules
#          tests/runtest.pl
#          tests/scf-opt/Makefile
#          tests/scf-opt2/Makefile
#          tests/scf-opt3/Makefile
#          tests/scf-opt4/Makefile
#          tests/scf-opt5/Makefile
#          tests/scf-opt6/Makefile
#          tests/scf-opt7/Makefile
#          tests/scf-opt9/Makefile
#          tests/scf-opt10/Makefile
#          tests/scf-opt11/Makefile
#          tests/scf-opt12/Makefile
#          tests/scf-opt-numer/Makefile
#          tests/scf-opt2-numer/Makefile
#          tests/scf-freq/Makefile
#          tests/scf-freq-first-3/Makefile
#          tests/scf-freq-first-5/Makefile
#          tests/scf-freq-none-3/Makefile
#          tests/scf-freq-none-5/Makefile
#          tests/scf-symm-fc-numer/Makefile
#          tests/scf-fc-numer/Makefile
#          tests/scf-polar/Makefile
#          tests/rhf-oeprop/Makefile
#          tests/rhf-stab/Makefile
#          tests/rhf-stab/rhf-stab.pl
#          tests/uhf-stab/Makefile
#          tests/uhf-stab/uhf-stab.pl
#          tests/rohf-stab/Makefile
#          tests/rohf-stab/rohf-stab.pl
#          tests/cis-sp/Makefile
#          tests/cisd-sp/Makefile
#          tests/cisd-opt-numer/Makefile
#          tests/fci-h2o/Makefile
#          tests/casscf-sp/Makefile
#          tests/casscf-sa-sp/Makefile
#          tests/casscf-fzc-sp/Makefile
#          tests/rasscf-sp/Makefile
#          tests/mp2-sp/Makefile
#          tests/mp2-direct-sp/Makefile
#          tests/mp2-oeprop/Makefile
#          tests/cc1/Makefile
#          tests/cc2/Makefile
#          tests/cc3/Makefile
#          tests/cc4/Makefile
#          tests/cc5/Makefile
#          tests/cc6/Makefile
#          tests/cc7/Makefile
#          tests/cc8/Makefile
#          tests/cc9/Makefile
#          tests/cc10/Makefile
#          tests/cc11/Makefile
#          tests/cc12/Makefile
#          tests/cc13/Makefile
#          tests/cc14/Makefile
#          tests/cc15/Makefile
#          tests/cc16/Makefile
#          tests/cc17/Makefile
#          tests/cc18/Makefile
#          tests/cc19/Makefile
#          tests/cc21/Makefile
#          tests/cc22/Makefile
#          tests/cc23/Makefile
#          tests/cc24/Makefile
#          tests/cc25/Makefile
#          tests/cc26/Makefile
#          tests/cc27/Makefile
#          tests/cc28/Makefile
#          tests/cc29/Makefile
#          tests/cc30/Makefile
#          tests/cc31/Makefile
#          tests/cc32/Makefile
#          tests/cc33/Makefile
#          tests/cc34/Makefile
#          tests/cc35/Makefile
#          tests/cc36/Makefile
#          tests/cc37/Makefile
#          tests/cc38/Makefile
#          tests/cc39/Makefile
#          tests/cc40/Makefile
#          tests/cc41/Makefile
#          tests/cc42/Makefile
#          tests/cc43/Makefile
#          tests/cc44/Makefile
#          tests/cc45/Makefile
#          tests/cc46/Makefile
#          tests/cc47/Makefile
#          tests/cc48/Makefile
#          tests/cc49/Makefile
#          tests/cc50/Makefile
#          tests/dboc-rhf1/Makefile
#          tests/dboc-rohf1/Makefile
#          tests/dboc-uhf1/Makefile
#          tests/dboc-rcisd1/Makefile
#          tests/dboc-rocisd1/Makefile
#          tests/extrema-zmat/Makefile
#          tests/extrema-zmat/extrema-zmat.pl
#          tests/extrema-deloc/Makefile
#          tests/extrema-deloc/extrema-deloc.pl
#          tests/mp2r12-sp1/Makefile
#          tests/mp2-opt-numer/Makefile
#          tests/mp2-opt/Makefile
#          tests/rhf-lindep1/Makefile
#          tests/rohf-lindep1/Makefile
#          tests/uhf-lindep1/Makefile
#          tests/psi_start/Makefile
#          tests/psi_start/psi_start.pl
#          src/bin/giao/Makefile

AC_OUTPUT
