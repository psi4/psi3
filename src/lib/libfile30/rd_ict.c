#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_ict():  Reads the transformation properties of the nuclei
**     under the operations allowed for the particular symmetry point group 
**     in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: int **ict  a matrix of integers. Each row corresponds 
**     to a particular symmetry operation, while each column corresponds to
**     a particular atom.  The value of ict[2][1], then, should be interpreted 
**     in the following manner: under the third symmetry operation of the 
**     relavant point group, the second atom is placed in the location
**     originally occupied by the atom with the index ict[2][1].
*/

int **file30_rd_ict(void)
{
  int i,nirreps,natom;
  PSI_FPTR ict_ptr;
  int **ict;

  nirreps = file30_rd_nirreps();
  natom = file30_rd_natom();
  ict_ptr = (PSI_FPTR) (info30_.mpoint[1] - 1)*sizeof(int);

  ict = (int **)malloc(sizeof(char *)*nirreps);
  for(i=0;i<nirreps;i++) 
    {
     ict[i] = (int *)malloc(sizeof(int)*natom);
     wreadw(info30_.filenum, (char *) ict[i], (int) sizeof(int)*natom,
	 ict_ptr, &ict_ptr);
    }

  return ict;
}
