#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void wwritw(int, char *, int, PSI_FPTR, PSI_FPTR *);

/* writes nlen bytes from array to file itape starting at current */
/* pointer location */

void swrit(itape,array,nlen)
   int itape, nlen;
   char *array;
   {
      PSI_FPTR start, end;

      start=ptr.wptr[itape];

      wwritw(itape,(char *) array,nlen,start,&end);
      ptr.wptr[itape]=((ptr.wptr[itape]-1+4096)/4096)*4096;
   }
