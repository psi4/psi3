/*!
  \file ict.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_ict():  Reads the transformation properties of the nuclei
**     under the operations allowed for the particular symmetry point group 
**     in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: ict = a matrix of integers. Each row corresponds 
**     to a particular symmetry operation, while each column corresponds to
**     a particular atom.  The value of ict[2][1], then, should be interpreted 
**     in the following manner: under the third symmetry operation of the 
**     relavant point group, the second atom is placed in the location
**     originally occupied by the atom with the index ict[2][1].
** \ingroup (CHKPT)
*/

int **chkpt_rd_ict(void)
{
  int i, natom, nirreps;
  int **ict;
  psio_address ptr;
  char *key;

  nirreps = chkpt_rd_nirreps();
  natom = chkpt_rd_natom();

  ptr = PSIO_ZERO;
  ict = (int **) malloc(sizeof(char *) * nirreps);
  key = chkpt_build_keyword("ICT Table");
  for(i=0; i < nirreps; i++) {
    ict[i] = (int *) malloc(sizeof(int) * natom);
    psio_read(PSIF_CHKPT, key, (char *) ict[i], natom*sizeof(int), ptr, &ptr);
  }
  free(key);

  return ict;
}


/*!
** chkpt_wt_ict():  Reads the transformation properties of the nuclei
**     under the operations allowed for the particular symmetry point group 
**     in which the molecule is considered.
**
**   arguments:
**   \param ict = a matrix of integers. Each row corresponds 
**     to a particular symmetry operation, while each column corresponds to
**     a particular atom.  The value of ict[2][1], then, should be interpreted 
**     in the following manner: under the third symmetry operation of the 
**     relavant point group, the second atom is placed in the location
**     originally occupied by the atom with the index ict[2][1].
**
**   returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_ict(int **ict)
{
  int i, natom, nirreps;
  psio_address ptr;
  char *key;

  nirreps = chkpt_rd_nirreps();
  natom = chkpt_rd_natom();

  key = chkpt_build_keyword("ICT Table");
  ptr = PSIO_ZERO;
  for(i=0; i < nirreps; i++)
    psio_write(PSIF_CHKPT, key, (char *) ict[i], natom*sizeof(int), 
               ptr, &ptr);
  free(key);
}
