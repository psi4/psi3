/*!
  \file atom_pos.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_atom_position() 
**  Reads in symmetry positions of atoms:
**    Possible values are as follows:
**	1   - atom in general position
**      2   - atom on c2z axis
**	4   - atom on c2y axis
**	8   - atom on c2x axis
**	16  - atom in the inversion center
**	32  - atom in the sigma_xy plane
**	64  - atom in the sigma_xz plane
**	128 - atom in the sigma_yz plane
**	This data is sufficient to define stabilizers of the nuclei.
**
**  returns: int *atom_position  an array of symmetry positions of atoms 
**  \ingroup (CHKPT)
*/

int *chkpt_rd_atom_position(void)
{
  int *atom_position, natom;
  char *key;

  natom = chkpt_rd_natom();
  atom_position = init_int_array(natom);

  key = chkpt_build_keyword("Atomic symm positions");
  psio_read_entry(PSIF_CHKPT, key, (char *) atom_position, natom*sizeof(int));
  free(key);

  return atom_position;
}


/*!
** chkpt_wt_atom_position() 
**
**  Writes out symmetry positions of atoms:
**    Possible values are as follows:
**	1   - atom in general position
**      2   - atom on c2z axis
**	4   - atom on c2y axis
**	8   - atom on c2x axis
**	16  - atom in the inversion center
**	32  - atom in the sigma_xy plane
**	64  - atom in the sigma_xz plane
**	128 - atom in the sigma_yz plane
**	This data is sufficient to define stabilizers of the nuclei.
**
**  \param atom_position = an array of symmetry positions of atoms
**
**  returns: none
**  \ingroup (CHKPT)
*/

void chkpt_wt_atom_position(int *atom_position)
{
  int natom;
  char *key;

  natom = chkpt_rd_natom();

  key = chkpt_build_keyword("Atomic symm positions");
  psio_write_entry(PSIF_CHKPT, key, (char *) atom_position, natom*sizeof(int));
  free(key);
}
