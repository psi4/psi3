/*!
** \file rfile.c
** \ingroup (CIOMR)
*/

#include "includes.h"
#include "pointers.h"

extern void ioopen_(int *);
extern void init_ptrs(void);

/*!
** rfile: open a binary file
**
** \param unit = file number
**
** \ingroup (CIOMR)
*/
void rfile(int unit)
{
       if (ptr.wptr == NULL) init_ptrs();

       sector = 1024;
       ptr.wptr[unit] = 0;
       ioopen_(&unit);
    }



