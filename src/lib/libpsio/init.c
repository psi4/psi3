/*!
** \file init.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "psio.h"

/* Definitions of global structs */
psio_ud *psio_unit;
psio_address PSIO_ZERO = {0,0};

#ifdef PSIO_STATS
ULI *psio_readlen;
ULI *psio_writlen;
#endif

/*!
** PSIO_INIT(): Allocates global memory needed by the I/O routines.
**
** No arguments.
*/

int psio_init(void)
{
  int i,j;

  psio_unit = (psio_ud *) malloc(sizeof(psio_ud)*PSIO_MAXUNIT);

#ifdef PSIO_STATS
  psio_readlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
  psio_writlen = (ULI *) malloc(sizeof(ULI) * PSIO_MAXUNIT);
#endif

  if(psio_unit == NULL) {
      fprintf(stderr, "Error in PSIO_INIT()!\n");
      exit(1);
    }

  for(i=0; i < PSIO_MAXUNIT; i++) {
#ifdef PSIO_STATS
      psio_readlen[i] = psio_writlen[i] = 0;
#endif      
      psio_unit[i].numvols = 0;
      for(j=0; j < PSIO_MAXVOL; j++) {
          psio_unit[i].vol[j].path = NULL;
          psio_unit[i].vol[j].stream = -1;
        }
      psio_unit[i].tocaddress.page = 0;
      psio_unit[i].tocaddress.offset = 0;
      psio_unit[i].toclen = 0;
      psio_unit[i].toc = NULL;
    }

  return(0);
}
