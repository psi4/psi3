
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:29:22  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

void int_pac(i,ib,j,jb,ipak,jpak)
   int *i,*ib,*j,*jb;
   unsigned int *ipak,*jpak;

   {
      *ipak = *ib;
      *ipak <<= 4;
      *ipak += *i;


      *jpak = *jb;
      *jpak <<= 4;
      *jpak += *j;
      }

void int_unpac(i,ib,j,jb,ipak,jpak)
   int *i,*ib,*j,*jb;
   unsigned int *ipak, *jpak;

   {
      *ib = *ipak >> 4;
      *i = *ipak & 15;

      *jb = *jpak >> 4;
      *j = *jpak & 15;
      }
