dnl
dnl  Local macros.  This are based on files provided with autoconf
dnl so the GNU GPL applies to them.
dnl
dnl
dnl See if the given fortran program compiles OK.  This requires the
dnl shell variable F77SUF to be defined.
dnl arg1: echo text
dnl arg2: optional code (main is already provided)
dnl arg3: additional compiler arguments
dnl arg4: success action
dnl arg5: fail action
dnl
define(AC_FC_COMPILE_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_MSG_CHECKING([for $1])]
)dnl
cat > conftest.$F77SUF <<EOF
      PROGRAM MAIN
      CALL X()
      STOP
      END

      SUBROUTINE X()
[$2]
      RETURN
      END
EOF
dnl Don't try to run the program, which would prevent cross-configuring.
if eval $FC [$3] -o conftest conftest.${F77SUF} ${FLIBS} > /dev/null  2>&1; then
  ifelse([$4], , :, [rm -rf conftest*
  $4
])
AC_MSG_RESULT(yes)
else
ifelse([$5], , , [rm -f conftest*
  $5
])dnl
AC_MSG_RESULT(no)
fi
rm -f conftest*]
)dnl
dnl
dnl See if the given C program is processed by the C compiler OK.  Does
dnl not necessarily compile.  If you want to compile give a -c to the compiler
dnl arguments.  To compile and link give a -o conftest.
dnl arg1: echo text
dnl arg2: optional code (main is already provided)
dnl arg3: additional compiler arguments
dnl arg4: success action
dnl arg5: fail action
dnl
define(AC_CC_PROCESS_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_MSG_CHECKING([for $1])]
)dnl
cat > conftest.c <<EOF
[$2]
main() {}
EOF
dnl Don't try to run the program, which would prevent cross-configuring.
if eval $CC [$3] conftest.c >/dev/null 2>&1; then
  ifelse([$4], , :, [rm -rf conftest*
  $4
])
ifelse([$5], , , [else
  rm -rf conftest*
  $5
])dnl
fi
rm -f conftest*
AC_MSG_RESULT(OK)]
)dnl
dnl
dnl  This checks for the existence of fortran libraries.  It is much
dnl like AC_HAVE_LIBRARY, except it uses AC_FC_COMPILE_CHECK
dnl 
define(AC_HAVE_FC_LIBRARY, [dnl
changequote(/,/)dnl
define(/AC_LIB_NAME/, dnl
patsubst(patsubst($1, /lib\([^\.]*\)\.a/, /\1/), /-l/, //))dnl
changequote([,])dnl
ac_save_FLIBS="${FLIBS}"
FLIBS="${FLIBS} -l[]AC_LIB_NAME[]"
ac_have_lib=""
AC_FC_COMPILE_CHECK([-l[]AC_LIB_NAME[]], , -o conftest, [ac_have_lib="1"])dnl
FLIBS="${ac_save_FLIBS}"
ifelse($#, 1, [dnl
if test -n "${ac_have_lib}"; then
   AC_DEFINE([HAVE_LIB]translit(AC_LIB_NAME, [a-z], [A-Z]))
   FLIBS="${FLIBS} -l[]AC_LIB_NAME[]"
fi
undefine(AC_LIB_NAME)dnl
], [dnl
if test -n "${ac_have_lib}"; then
   :; $2
else
   :; $3
fi
])])dnl
dnl
dnl try to determine what main is called, and how to link
dnl c and fortran codes
dnl shell variable F77SUF to be defined.
dnl arg1: echo text
dnl arg2: additional compiler arguments
dnl
define(AC_FC_LINKAGE_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_MSG_CHECKING([for $1])]
)dnl
cat > conftest.$F77SUF <<EOF
      PROGRAM FOO
      STOP
      END
      SUBROUTINE FOO2
      RETURN
      END
EOF
dnl check to see if this is gnu egrep which uses -q for silent mode,
dnl or a more generic one which uses -s
if eval egrep -q FOO2 conftest.$F77SUF > /dev/null 2>&1; then
  ac_egrep_silent=-q
else
  ac_egrep_silent=-s
fi
dnl compile conftest.f and look for main and foo2 symbols
$FC [$2] -c conftest.${F77SUF} ${FLIBS} > /dev/null  2>&1
if eval nm -p conftest.o | egrep $ac_egrep_silent MAIN__; then
  MAIN_FUNC=MAIN__
elif eval nm -p conftest.o | egrep $ac_egrep_silent MAIN_; then
  MAIN_FUNC=MAIN_
elif eval nm -p conftest.o | egrep $ac_egrep_silent MAIN; then
  MAIN_FUNC=MAIN
elif eval nm -p conftest.o | egrep $ac_egrep_silent main__; then
  MAIN_FUNC=main__
elif eval nm -p conftest.o | egrep $ac_egrep_silent main_; then
  MAIN_FUNC=main_
elif eval nm -p conftest.o | egrep $ac_egrep_silent main; then
  MAIN_FUNC=main
fi
AC_MSG_RESULT($1 is $MAIN_FUNC)
CDEF="$CDEF -DMAIN_FUNC=$MAIN_FUNC"
dnl
if eval nm -p conftest.o | egrep $ac_egrep_silent foo2_; then
  CDEF="$CDEF -DFCLINK=1"
  CXXDEF="$CXXDEF -DFCLINK=1"
elif eval nm -p conftest.o | egrep $ac_egrep_silent foo2; then
  CDEF="$CDEF -DFCLINK=2"
  CXXDEF="$CXXDEF -DFCLINK=2"
elif eval nm -p conftest.o | egrep $ac_egrep_silent FOO2; then
  CDEF="$CDEF -DFCLINK=3"
  CXXDEF="$CXXDEF -DFCLINK=3"
fi
rm -f conftest*]
)dnl
dnl
dnl try to determine the sizes of some c basic types
dnl arg1: echo text
dnl
define(AC_C_SIZES_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_CHECKING([for $1])]
)dnl
cat > conftest.c <<EOF
#include <stdio.h>
main()
{
  printf("#!/bin/sh\n");
  printf("ac_sizeof_int=%d\n",sizeof(int));
  printf("ac_sizeof_long=%d\n",sizeof(long));
  printf("ac_sizeof_long_long=%d\n",sizeof(long long));
  printf("ac_sizeof_double=%d\n",sizeof(double));
  printf("ac_sizeof_long_double=%d\n",sizeof(long double));
  printf("ac_sizeof_pointer=%d\n",sizeof(int*));
  return 0;
}
EOF
$CC conftest.c -o conftest
./conftest > conftest.sh
. conftest.sh
cat > conftest.$F77SUF <<EOF
      program foo
      real*8 a(2)
      double precision d(2)
      integer i(2)
      integer aa1 , aa2, ai1, ai2, ad1, ad2
      aa1 = %loc(a(1))
      aa2 = %loc(a(2))
      ai1 = %loc(i(1))
      ai2 = %loc(i(2))
      ad1 = %loc(d(1))
      ad2 = %loc(d(2))
      write(6,12)
      write(6,10) aa2-aa1
      write(6,11) ai2-ai1
      write(6,13) ad2-ad1
  10  format('ac_sizeof_real=',i1)
  11  format('ac_sizeof_integer=',i1)
  12  format('#!/bin/sh')
  13  format('ac_sizeof_dprec=',i1)
      stop
      end
EOF
if eval $FC $FOTH conftest.${F77SUF} -o conftest > /dev/null 2>&1; then
  ./conftest > conftest.sh
  . conftest.sh
else
  ac_sizeof_real=$ac_sizeof_double
  ac_sizeof_dprec=$ac_sizeof_double
  ac_sizeof_integer=$ac_sizeof_int
fi
rm -f conftest*]
)dnl
dnl
dnl Check for 64-bit capable IBM RS/6000
dnl arg1 : echo text
dnl
define(AC_RS64_CHECK,
[AC_PROVIDE([$0])dnl
ifelse([$1], , , [AC_MSG_CHECKING([for $1])]
)dnl
cat > conftest.c <<EOF
#include <stdio.h>
#include "/usr/include/sys/systemcfg.h"

int main()
{
  printf("#!/bin/sh\n");
  if (__power_64())
    printf("bitwidth64=yes\n");
  else
    printf("bitwidth64=no\n");
  exit(0);
}
EOF
$CC conftest.c -o conftest
./conftest > conftest.sh
. conftest.sh
rm -f conftest*
AC_MSG_RESULT("$bitwidth64")]
)dnl

dnl
dnl Check for LaTeX
dnl
define(AC_PROG_LATEX,
[AC_PROVIDE([$0])
AC_CHECK_PROG(LATEX, latex, latex)]
)dnl

dnl
dnl Check for dvips
dnl
define(AC_PROG_DVIPS,
[AC_PROVIDE([$0])
AC_CHECK_PROG(DVIPS, dvips, dvips)]
)dnl


dnl
dnl Check for LaTeX2HTML
dnl
define(AC_PROG_LATEX2HTML,
[AC_PROVIDE([$0])
AC_CHECK_PROG(LATEX2HTML, latex2html, latex2html)]
)dnl

