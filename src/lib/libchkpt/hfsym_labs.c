/*!
  \file hfsym_labs.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_hfsym_labs(): Read in the symmetry labels for all irreps in the 
**   point group in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: char **hfsym_labs   an array of labels (strings) which denote
**      the irreps which have basis functions (in Cotton ordering).  For DZ or
**      STO water, for example, in C2v symmetry, this would be an array of 
**      three labels: "A1", "B1", and "B2".
*/

char **chkpt_rd_hfsym_labs(void)
{
  int i, nsymhf, nirreps;
  char **hfsym_labs;
  psio_address ptr;

  nsymhf = chkpt_rd_nsymhf();
  nirreps = chkpt_rd_nirreps();

  ptr = PSIO_ZERO;
  hfsym_labs = (char **)malloc(sizeof(char *)*nirreps);
  for(i=0;i<nsymhf;i++) {
    hfsym_labs[i] = (char *) malloc(4*sizeof(char));
    psio_read(PSIF_CHKPT, "::HF irrep labels", (char *) hfsym_labs[i], 
	      4*sizeof(char), ptr, &ptr);
  }

  return hfsym_labs;
}
