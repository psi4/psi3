#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

double d1diag(void)
{
  int h, nirreps, i;
  double **T, **C, *E, max;
  struct oe_dpdfile T1;

  nirreps = moinfo.nirreps;
  max = 0.0;

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  for(h=0; h < nirreps; h++) {
      if(T1.params->rowtot[h]) {
         T = block_matrix(T1.params->rowtot[h], T1.params->rowtot[h]);

         newmm(T1.matrix[h], 0, T1.matrix[h], 1, T, T1.params->rowtot[h],
   	       T1.params->coltot[h], T1.params->rowtot[h], 1.0, 0.0);

         E = init_array(T1.params->rowtot[h]);
         C = block_matrix(T1.params->rowtot[h], T1.params->rowtot[h]);
         sq_rsp(T1.params->rowtot[h], T1.params->rowtot[h], T, E, 0, C, 1e-12);

         /* Find maximum eigenvalue of T */
         for(i=0; i < T1.params->rowtot[h]; i++) if(E[i] > max) max = E[i];
	     
         free_block(T);
         free_block(C);
         free(E);
       }
    }

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);

  max = sqrt(max);

  return max;
}
