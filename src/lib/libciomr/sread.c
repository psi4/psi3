#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void wreadw(int, char *, int, PSI_FPTR, PSI_FPTR *);

/* reads nlen bytes from itape into array starting at current pointer */
/* location */

void sread(itape,array,nlen)
   int itape, nlen;
   char *array;
   {
      PSI_FPTR start, end;

      start=ptr.wptr[itape];

      wreadw(itape,(char *) array,nlen,start,&end);
      ptr.wptr[itape]=((ptr.wptr[itape]-1+4096)/4096)*4096;
   }
