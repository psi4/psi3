#include <math.h>
#include <libdpd/dpd.h>
#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"
#define EXTERN
#include "globals.h"

/* onepdm(): Computes the one-particle density matrix for MP2 wave functions.
**
** The spin-orbital expression for the OPDM:
**
** docc-docc   D_ij = -1/2 T_ik^ab T_jk^ab
**
** virt-virt   D_ab =  1/2 T_ij^ac T_ij^bc
**
*/

void opdm(void)
{
  double trace=0.0;
  dpdfile2 D;
  dpdbuf4 T2A, T2B;
  int h, i, j, a, b;
  int I, J, A, B;
  double **OPDM;
  
  /*
  ** D_IJ = -2 * SUM(k,a,b) T_ik^ab * T_jk^ab + SUM(k,a,b) T_ik^ab * T_jk^ba 
  */
  
  /* Form tIjAb and tIjBa */
  dpd_buf4_init(&T2A, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_scmcopy(&T2A, CC_TMP0, "2 tIjAb - tIjBa", 2);
  dpd_buf4_sort_axpy(&T2A, CC_TMP0, pqsr, 0, 5, "2 tIjAb - tIjBa", -1);
  dpd_buf4_close(&T2A);
  
  /* DOCC-DOCC Block of OPDM(2) */
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&T2A, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TMP0, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_contract442(&T2A, &T2B, &D, 0, 0, -2, 0);
  dpd_buf4_close(&T2B);
  dpd_buf4_close(&T2A);
  trace += dpd_file2_trace(&D);
  dpd_file2_close(&D);
 
  /*
  ** D_AB = 2 * SUM(c,i,j) T_ij^ac T_ij^bc - SUM(c,i,j) T_ij^ac T_ji^bc 
  */

  /* VIRT-VIRT Block of OPDM(2) */
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&T2A, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TMP0, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_contract442(&T2A, &T2B, &D, 3, 3, 2, 0);
  dpd_buf4_close(&T2B);
  dpd_buf4_close(&T2A);
  trace += dpd_file2_trace(&D);
  dpd_file2_close(&D);
  
  fprintf(outfile, "\n\tTrace of OPDM(2) \t=\t  %.10f\n", fabs(trace));

  OPDM = block_matrix(mo.nmo,mo.nmo);  

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  
  trace = 0.0;
  
  for(h=0; h<mo.nirreps; h++) {
    for(i=0; i<mo.doccpi[h]; i++) {
        I = mo.qt_docc[mo.docc_off[h] + i];
      for(j=0; j<mo.doccpi[h]; j++) {
	J = mo.qt_docc[mo.docc_off[h] + j];
	OPDM[I][J] += D.matrix[h][i][j];
	/* Add OPDM(HF) to OPDM(2) */
	if(I==J) {
	  OPDM[I][J] += 2;
	  trace += OPDM[I][J];
	}
      }
    }
  }
   
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
    
  for(h=0; h<mo.nirreps; h++) {
    for(a=0; a<mo.virtpi[h]; a++) {
        A = mo.qt_virt[mo.virt_off[h] + a];
      for(b=0; b<mo.virtpi[h]; b++) {
	B = mo.qt_virt[mo.virt_off[h] + b];
	OPDM[A][B] += D.matrix[h][a][b];
	if(A==B) {
	  trace += OPDM[A][B];
	}
      }
    }
  }

  dpd_file2_mat_close(&D);
  fprintf(outfile, "\tTrace of MP2 OPDM \t=\t %.10f\n", trace);
  dpd_file2_close(&D);

  if (params.opdm_print || params.print > 5) {
    fprintf(outfile, "\n\tMP2 OPDM (MO)\n");
    print_mat(OPDM,mo.nmo,mo.nmo,outfile);
  }

  if (params.opdm_write) {
    psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_write_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) OPDM[0],
                     sizeof(double)*mo.nmo*mo.nmo);
    psio_close(PSIF_MO_OPDM, 1);
  }
  
  free_block(OPDM);
}
