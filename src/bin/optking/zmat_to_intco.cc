/*! \file 
    \ingroup (OPTKING)
    \brief ZMAT_TO_INTCO() determine simples from z-matrix coordinates 
  - may not work for dummy atoms
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>

#define EXTERN
#include "opt.h"
#undef EXTERN

namespace psi { namespace optking {

void zmat_to_intco() {
  int i, first, a, b, c, d, cnt = 0, natom;
  char **felement;
  char buf[2];
  struct z_entry *zmat;

  natom = optinfo.natom;
  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  chkpt_close();

/*
  for (i=0;i<natom;++i) {
    fprintf(outfile,"%d %d %d %d\n",  i, zmat[i].bond_atom, zmat[i].angle_atom, zmat[i].tors_atom);
    fprintf(outfile,"%20.10lf %20.10lf %20.10lf\n", zmat[i].bond_val, zmat[i].angle_val, zmat[i].tors_val);
    fprintf(outfile,"%s %s %s \n", zmat[i].bond_label, zmat[i].angle_label, zmat[i].tors_label);
   }
*/

  ffile(&fp_intco,"intco.dat",0);
  fprintf(fp_intco,"intco: (\n");

  fprintf(fp_intco,"  stre = (\n");
  for (i=1; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    swap(&a, &b);
    fprintf(fp_intco, "    (%d %d %d)\n",++cnt, a, b);
  }
  fprintf(fp_intco,"  )\n");
  
  fprintf(fp_intco,"  bend = (\n");
  for (i=2; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    swap(&a, &c);
    if (zmat[i].angle_val != 180.0) {
      fprintf(fp_intco, "    (%d %d %d %d)\n",++cnt, a, b, c);
    }
  }
  fprintf(fp_intco,"  )\n");

  fprintf(fp_intco,"  tors = (\n");
  for (i=3; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    d = zmat[i].tors_atom;
    swap_tors(&a, &b, &c, &d);
    if (zmat[i].angle_val != 180.0)
      fprintf(fp_intco, "    (%d %d %d %d %d)\n",++cnt, a, b, c, d);
  }
  fprintf(fp_intco,"  )\n");

  first = 1;
  for (i=2; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    swap(&a, &c);
    if (zmat[i].angle_val == 180.0) {
      if (first) {
        fprintf(fp_intco,"  lin1 = (\n");
        first = 0;
      }
      fprintf(fp_intco, "    (%d %d %d %d)\n",++cnt, a, b, c);
    }
  }
  if (!first) fprintf(fp_intco,"  )\n");

  first = 1;
  for (i=2; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    swap(&a, &c);
    if (zmat[i].angle_val == 180.0) {
      if (first) {
        fprintf(fp_intco,"  lin2 = (\n");
        first = 0;
      }
      fprintf(fp_intco, "    (%d %d %d %d)\n",++cnt, a, b, c);
    }
  }
  if (!first) fprintf(fp_intco,"  )\n");
  fprintf(fp_intco,")");
  fclose(fp_intco);

  ffile(&fp_intco,"fintco.dat",0);

  /* write out coordinates to be frozen */
  fprintf(fp_intco,"fixed_intco: (\n");
  cnt = 0;

  fprintf(fp_intco,"  stre = (\n");
  for (i=1; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    swap(&a, &b);
    if (zmat[i].bond_opt == 0)
      fprintf(fp_intco, "    ( %d %d)\n", a, b);
  }
  fprintf(fp_intco,"  )\n");

  fprintf(fp_intco,"  bend = (\n");
  for (i=2; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    swap(&a, &c);
    if (zmat[i].angle_val != 180.0) {
      if (zmat[i].angle_opt == 0)
        fprintf(fp_intco, "    ( %d %d %d)\n", a, b, c);
    }
  }
  fprintf(fp_intco,"  )\n");

  fprintf(fp_intco,"  tors = (\n");
  for (i=3; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    d = zmat[i].tors_atom;
    swap_tors(&a, &b, &c, &d);
    if (zmat[i].angle_val != 180.0) {
      if (zmat[i].tors_opt == 0)
        fprintf(fp_intco, "    ( %d %d %d %d)\n", a, b, c, d);
    }
  }
  fprintf(fp_intco,"  )\n");

  first = 1;
  for (i=2; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    swap(&a, &c);
    if (zmat[i].angle_val == 180.0) {
      if (zmat[i].angle_opt == 0) {
        if (first) {
          fprintf(fp_intco,"  lin1 = (\n");
          first = 0;
        }
        fprintf(fp_intco, "    ( %d %d %d)\n", a, b, c);
      }
    }
  }
  if (!first) fprintf(fp_intco,"  )\n");

  first = 1;
  for (i=2; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    c = zmat[i].angle_atom;
    swap(&a, &c);
    if (zmat[i].angle_val == 180.0) {
      if (zmat[i].angle_opt == 0) {
        if (first) {
          fprintf(fp_intco,"  lin2 = (\n");
          first = 0;
        }
        fprintf(fp_intco, "    ( %d %d %d)\n", a, b, c);
      }
    }
  }
  if (!first) fprintf(fp_intco,"  )\n");

  if (cnt > 0)
    optinfo.constraints_present = 1;

  fprintf(fp_intco,")");
  fclose(fp_intco);
  return;
}

}} /* namespace psi::optking */
