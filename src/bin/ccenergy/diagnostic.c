#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double diagnostic(void)
{
  int h, nirreps, Gi, Ga;
  int i, a, I, A, row, col;
  int open_shell, num_elec;
  int *occpi, *virtpi;
  int *occ_sym, *vir_sym;
  int *clsdpi, *uoccpi;
  int *openpi;
  double t1diag;
  struct oe_dpdfile T1A, T1B;

  nirreps = moinfo.nirreps;
  clsdpi = moinfo.clsdpi; uoccpi = moinfo.uoccpi;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Check for open shells */
  open_shell = 0;
  for(h=0; h < nirreps; h++)
      if(openpi[h]) open_shell = 1;

  /* Compute the number of electrons */
  for(h=0,num_elec=0; h < nirreps; h++)
      num_elec += (2 * clsdpi[h] + openpi[h]);

  dpd_oe_file_init(&T1A, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1A);
  dpd_oe_file_mat_rd(&T1A, 0, outfile);
  dpd_oe_file_init(&T1B, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1B);
  dpd_oe_file_mat_rd(&T1B, 0, outfile);

  /* Closed-shell diagnostic --- T1A should be equal to T1B */
  t1diag = 0.0;
  if(!open_shell) {
      for(h=0; h < nirreps; h++) {

	  for(row=0; row < T1A.params->rowtot[h]; row++) 
	      for(col=0; col < T1A.params->coltot[h]; col++) 
		  t1diag += (T1A.matrix[h][row][col] * T1A.matrix[h][row][col]);
	}


      t1diag /= num_elec;
      t1diag = sqrt(t1diag);
    }
  /* Open-shell diagnostic */
  else {
      for(h=0; h < nirreps; h++) {

	  /* Loop over docc-docc components */
	  for(i=0; i < (occpi[h] - openpi[h]); i++) {
	      for(a=0; a < (virtpi[h] - openpi[h]); a++) {

		  t1diag += (T1A.matrix[h][i][a] + T1B.matrix[h][i][a]) *
			    (T1A.matrix[h][i][a] + T1B.matrix[h][i][a]);
		}
	    }

	  /* Loop over docc-socc components */
	  for(i=0; i < (occpi[h] - openpi[h]); i++) {
	      for(a=0; a < openpi[h]; a++) {
		  A = a + uoccpi[h];

		  t1diag += 2 * T1B.matrix[h][i][A] * T1B.matrix[h][i][A];
		}
	    }

	  /* Loop over socc-docc components */
	  for(i=0; i < openpi[h]; i++) {
	      I = i + clsdpi[h];
	      for(a=0; a < (virtpi[h] - openpi[h]); a++) {
		  
		  t1diag += 2 * T1A.matrix[h][I][a] * T1A.matrix[h][I][a];
		}
	    }
	}

      t1diag /= num_elec;
      t1diag = sqrt(t1diag);
      t1diag *= 0.5;
    }
	      
  dpd_oe_file_mat_close(&T1A);
  dpd_oe_file_close(&T1A);
  dpd_oe_file_mat_close(&T1B);
  dpd_oe_file_close(&T1B);

  return t1diag;
}
