#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void rhf_twopdm(void);
void uhf_twopdm(void);

void twopdm(void)
{
  if(params.ref == 0) rhf_twopdm();
  else if(params.ref == 2) uhf_twopdm();
}

void rhf_twopdm(void)
{
  int h,i,j,a,b,I,J,A,B,IJ,AB,Gi,Gj,Ga,II,AA;
  double trace = 0.0;
  dpdbuf4 T;
  dpdbuf4 G;

  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_copy(&T, CC_GAMMA, "GIjAb");
  dpd_buf4_close(&T);

  /*
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < mo.nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G,h);
    dpd_buf4_mat_irrep_rd(&G,h);
    for(i=0; i < mo.occpi[h]; i++) {
      I = mo.occ_off[h] + i;
      II = G.params->rowidx[I][I];
      for(a=0; a < mo.virpi[h]; a++) {
        A = mo.vir_off[h] + a;
        AA = G.params->colidx[A][A];

        trace += G.matrix[h][II][AA];
      }
    }
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  fprintf(outfile,"\n\tTrace of TPDM(2)        = %20.15f\n", fabs(trace));
  */
}

void uhf_twopdm(void)
{

}
