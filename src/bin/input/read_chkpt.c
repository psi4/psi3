#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include "input.h"
#include "global.h"
#include "defines.h"

void read_chkpt_geom()
{
  int i, errcod;
  int atom;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label;

  chkpt_init(PSIO_OPEN_OLD);
  num_atoms = chkpt_rd_natom();
  num_allatoms = chkpt_rd_nallatom();
  if (num_atoms == 0)
    punt("GEOMETRY in the checkpoint file is empty!");
  else if (num_atoms > MAXATOM)
    punt("There are more atoms than allowed!");

  /* Grab subgroup and get rid of a possible blank */
  subgroup = chkpt_rd_sym_label();
  if (subgroup[2] == ' ') subgroup[2] = '\0';

  full_geom = chkpt_rd_fgeom();
  geometry = (double **) malloc(num_atoms*sizeof(double *));
  atom_dummy = chkpt_rd_atom_dummy();
  atom = 0;
  for(i=0;i<num_allatoms;i++)
    if (!atom_dummy[i]) {
      geometry[atom] = full_geom[i];
      ++atom;
    }
  nuclear_charges = chkpt_rd_zvals();
  full_element = chkpt_rd_felement();
  chkpt_close();

  element = (char **) malloc(sizeof(char *)*num_atoms);
  elemsymb_charges = init_array(num_atoms);
  for(i=0;i<num_atoms;i++) {
    element[i] = elem_name[(int)nuclear_charges[i]];
    elemsymb_charges[i] = nuclear_charges[i];
  }

  return;
}


