/*** ZMAT_TO_INTCO() determine simples and salcs (someday?) from z-matrix ***/ 

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

/*** ZMAT_TO_INTCO if simples are not already there and zmat_simples
 * is not turned off, then generate simple internals from zmatrix ***/
void zmat_to_intco() {
  int i, a, b, c, d, cnt = 0;
  int nallatom, natom;
  char **felement;
  char buf[2];
  struct z_entry *zmat;

  nallatom = optinfo.nallatom;
  natom = optinfo.natom;

  if ( !(optinfo.simples_present) && (optinfo.zmat_simples) ) {
    chkpt_init(PSIO_OPEN_OLD);
    zmat = chkpt_rd_zmat();
    chkpt_close();

    ffile(&fp_intco,"intco.dat",0);
    fprintf(fp_intco,"intco: (\n");

    fprintf(fp_intco,"  stre = (\n");
    for (i=1; i<nallatom; ++i) {
      a = i+1;
      b = zmat[i].bond_atom;
      swap(&a, &b);
      if (zmat[i].bond_opt)
        fprintf(fp_intco, "    (%d %d %d)\n",++cnt, a, b);
    }
    fprintf(fp_intco,"  )\n");
  
    fprintf(fp_intco,"  bend = (\n");
    for (i=2; i<nallatom; ++i) {
      a = i+1;
      b = zmat[i].bond_atom;
      c = zmat[i].angle_atom;
      swap(&a, &c);
      if (zmat[i].angle_opt)
        fprintf(fp_intco, "    (%d %d %d %d)\n",++cnt, a, b, c);
    }
    fprintf(fp_intco,"  )\n");
  
    fprintf(fp_intco,"  tors = (\n");
    for (i=3; i<nallatom; ++i) {
      a = i+1;
      b = zmat[i].bond_atom;
      c = zmat[i].angle_atom;
      d = zmat[i].tors_atom;
      swap_tors(&a, &b, &c, &d);
      if (zmat[i].tors_opt)
        fprintf(fp_intco, "    (%d %d %d %d %d)\n",++cnt, a, b, c, d);
    }
    fprintf(fp_intco,"  )\n");
  
    fprintf(fp_intco,")");
    fclose(fp_intco);
  }

  // free zmat?

  return;
}

