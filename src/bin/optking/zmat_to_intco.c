/*** ZMAT_TO_INTCO() determine simples from z-matrix coordinates 
  - not intended to work for dummy atoms ***/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>

#define EXTERN
#define C_EXTERN
#define C_CODE
#include "opt.h"
#undef C_CODE
#undef C_EXTERN
#undef EXTERN

void zmat_to_intco() {
  int i, first, a, b, c, d, cnt = 0, natom;
  char **felement;
  char buf[2];
  struct z_entry *zmat;

  natom = optinfo.natom;
  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  chkpt_close();

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

  /* write out coordinates to be frozen */
  ffile(&fp_intco,"fintco.dat",0);
  fprintf(fp_intco,"intco: (\n");
  cnt = 0;

  fprintf(fp_intco,"  stre = (\n");
  for (i=1; i<natom; ++i) {
    a = i+1;
    b = zmat[i].bond_atom;
    swap(&a, &b);
    if (zmat[i].bond_opt == 0)
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
      if (zmat[i].angle_opt == 0)
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
    if (zmat[i].angle_val != 180.0) {
      if (zmat[i].tors_opt == 0)
        fprintf(fp_intco, "    (%d %d %d %d %d)\n",++cnt, a, b, c, d);
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
        fprintf(fp_intco, "    (%d %d %d %d)\n",++cnt, a, b, c);
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
        fprintf(fp_intco, "    (%d %d %d %d)\n",++cnt, a, b, c);
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
