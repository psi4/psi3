#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern PSI_FPTR sec2i(int);
extern void wwritw(int, char *, int, PSI_FPTR, PSI_FPTR *);

/* writes nlen bytes from array to file itape starting at sector irec */

void rwrit(itape,array,nlen,irec)
   int itape, nlen, irec;
   char *array;
   {
    PSI_FPTR ipos, junk;

    ipos = sec2i(--irec);
    wwritw(itape,(char *) array,nlen,ipos,&junk);
   }
