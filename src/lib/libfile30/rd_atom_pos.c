/*!
  \file rd_atom_pos.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_atom_position(): Reads in symmetry positions of atoms:
**	Possible values are as follows:
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
**  takes no arguments.
**
**  returns: int *atom_position  an array of symmetry positions of atoms 
*/


int *file30_rd_atom_position(void)
{
  int *atom_position;
  int natom;
  PSI_FPTR junk;
  PSI_FPTR ap_ptr;

  natom = file30_rd_natom();
  ap_ptr = (PSI_FPTR) (info30_.mpoint[38] - 1)*sizeof(int);

  atom_position = init_int_array(natom);

  wreadw(info30_.filenum, (char *) atom_position, (int) natom*sizeof(int),
	 ap_ptr, &junk);

  return atom_position;
}
