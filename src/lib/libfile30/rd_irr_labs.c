#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_irr_labs(): Read in the symmetry labels for all irreps in the 
**   point group in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: char **irr_labs   an array of labels (strings) which denote
**      the irreps for the point group	in which the molecule is considered,
**      _regardless_ of whether there exist any symmetry orbitals which 
**      transform as that irrep.  
*/

char **file30_rd_irr_labs(void)
{
  int i,nirreps;
  PSI_FPTR irr_labs_ptr;
  char **irr_labs;

  nirreps = file30_rd_nirreps();
  irr_labs_ptr = (PSI_FPTR) (info30_.mpoint[15] - 1)*sizeof(int);

  irr_labs = (char **)malloc(sizeof(char *)*nirreps);
  for(i=0;i<nirreps;i++) 
    {
     irr_labs[i] = (char *)malloc(sizeof(int));
     wreadw(info30_.filenum, (char *) irr_labs[i], (int) sizeof(int),
	 irr_labs_ptr, &irr_labs_ptr);
     irr_labs[i][3] = '\0';
    }

  return irr_labs;
}
