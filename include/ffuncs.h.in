
#ifndef _ffuncs_h
#define _ffuncs_h

#include <psiconfig.h>

#if FCLINK==1
# define C_QABORT qabort_
# define C_MABORT mabort_
#elif FCLINK==2
# define C_QABORT qabort
# define C_MABORT mabort
#else
# define C_QABORT QABORT
# define C_MABORT MABORT
#endif

extern void C_QABORT(void);
extern void C_MABORT(void);

#endif /* _ffuncs_h */
