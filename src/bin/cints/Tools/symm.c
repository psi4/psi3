#include<stdio.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<file30.h>

#include<libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static void init_dp_table(void);

void init_symmetry()
{
  int i, j, count;
  
  Symmetry.symlabel = file30_rd_sym_label();
  Symmetry.nirreps = file30_rd_nirreps();
  Symmetry.num_so = file30_rd_nso();
  Symmetry.num_unique_shells = file30_rd_num_unique_shell();

  Symmetry.atom_positions = file30_rd_atom_position();
  Symmetry.us2s = file30_rd_us2s();
  Symmetry.sopi = file30_rd_sopi();
  Symmetry.sym_oper = file30_rd_symoper();
  Symmetry.irr_labels = file30_rd_irr_labs();

  if (Symmetry.nirreps) {
  /* Symmetry.dp_table = */ init_dp_table();
#if SCF_ONLY
    if (Symmetry.nirreps < 4)
      UserOptions.scf_only = 0;
    if (UserOptions.scf_only) {
      Symmetry.so2symblk = init_int_array(Symmetry.num_so);
      count = 0;
      for(i=0;i<Symmetry.nirreps;i++)
	for(j=0;j<Symmetry.sopi[i];j++)
	  Symmetry.so2symblk[count++] = i;
    }
#endif
  }

  return;
}


void cleanup_symmetry()
{
  if (Symmetry.nirreps > 1)
    free_int_matrix(Symmetry.dp_table,Symmetry.nirreps);
  free(Symmetry.sopi);
  free_block(Symmetry.usotao);
  free(Symmetry.us2s);
  free(Symmetry.atom_positions);

  return;
}


/*----------------------------------------------------------------------
  Compute direct product multiplication table for the given point group
  NOTE: This matrix is really pointless at the moment, I left it here
  just in case
 ----------------------------------------------------------------------*/
void init_dp_table(void)
{
  int i,j;

  Symmetry.dp_table = init_int_matrix(Symmetry.nirreps,
				      Symmetry.nirreps);
  for(i=0;i<Symmetry.nirreps;i++)
    for(j=0;j<i;j++) {
      /*------------------------------------------------------------
	The line below works only in a case of Abelian point group!
	Have to do honest multiplication of rows of character table
	if non-Abelian groups to be used
       ------------------------------------------------------------*/
      Symmetry.dp_table[i][j] = Symmetry.dp_table[j][i] = i ^ j;
    }

  return;
}
