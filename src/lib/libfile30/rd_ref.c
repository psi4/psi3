#include "file30.h"
#include "file30.gbl"

/*
** rd_ref():  Reads the reference type from the flag in file30.
**   0 = RHF | 1 = UHF | 2 = ROHF | 3 = TCSCF 
**
**   takes no arguments.
**
**   returns: int refnum   number indicating the reference.
*/

int file30_rd_ref(void)
{
    return info30_.mconst[46];
}
