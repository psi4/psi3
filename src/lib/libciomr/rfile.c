/*!
** \file rfile.c
** \ingroup (CIOMR)
*/

/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/19 21:48:06  sherrill
/* Remove some unused functions and do doxygen markup of libciomr.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:22  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.2  1999/11/01 20:10:58  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.1  1991/06/15 18:29:49  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

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



