/*** PRINT_ZMAT() printout z-matrix ***/ 

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
void print_zmat(FILE *outfile) {
  int i, a, b, c, d, cnt = 0;
  int nallatom, natom;
  char **felement;
  char buf[2];
  double *zvals;
  struct z_entry *zmat;

  nallatom = optinfo.nallatom;
  natom = optinfo.natom;

  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  felement = chkpt_rd_felement();
  chkpt_close();

  for (i=0; i<nallatom; ++i) {
    if (i == 0) {
      fprintf(outfile,"    ( %s )\n", felements[i]);
    }
    else if (i == 1) {
      fprintf(outfile,"    ( %s", felements[i]);
      fprintf(outfile,")\n");
    }
    else if (i == 2) {
      fprintf(outfile,"    ( %s", felements[i]);
      fprintf(outfile,")\n");
    }
    else if (i >=3) {
      fprintf(outfile,"    ( %s", felements[i]);
      fprintf(outfile,")\n");
    }
  }

      a = i+1;
      b = zmat[i].bond_atom;
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


  return;
}

