#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void wreadw(int, char *, int, PSI_FPTR, PSI_FPTR *);
extern PSI_FPTR sec2i(int);

/* reads nlen bytes into array starting at irec */

void rread(itape,array,nlen,irec)
   int itape, nlen, irec;
   char *array;
   {
      PSI_FPTR ipos, junk;

      ipos = sec2i(--irec);
      wreadw(itape,(char *) array,nlen,ipos,&junk);
   }
