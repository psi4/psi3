#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void rhf_mohess(void);
void uhf_mohess(void);

void mohess(void)
{
  if(params.ref == 0) rhf_mohess();
  else if(params.ref == 2) uhf_mohess();
}

void rhf_mohess(void)
{
  dpdbuf4 Amat;  /* MO Hessian */
  dpdbuf4 D;     /* Two-electron integral */
  dpdbuf4 C;     /* Two-electron integral */
  dpdfile2 fIJ;  /* Occ-Occ Fock matrix */
  dpdfile2 fAB;  /* Vir-Vir Fock matrix */
  int h, nirreps;
  int i, j, a, b, ai, bj;
  int I, J, A, B;
  int Asym, Bsym, Isym, Jsym;
  
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_MISC, rpsq, 11, 11, "A(AI,BJ)");
  dpd_buf4_close(&D);

  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_scm(&Amat, 4.0);
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort_axpy(&D, CC_MISC, sprq, 11, 11, "A(AI,BJ)", -1.0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_sort_axpy(&C, CC_MISC, qpsr, 11, 11, "A(AI,BJ)", -1.0);
  dpd_buf4_close(&C);
  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  dpd_buf4_close(&Amat);

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);

  nirreps = mo.nirreps;

  dpd_buf4_init(&Amat, CC_MISC, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);
    for(ai=0; ai < Amat.params->rowtot[h]; ai++) {
      a = Amat.params->roworb[h][ai][0];
      i = Amat.params->roworb[h][ai][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(bj=0; bj < Amat.params->coltot[h]; bj++) {
        b = Amat.params->colorb[h][bj][0];
        j = Amat.params->colorb[h][bj][1];
        B = fAB.params->colidx[b];
        J = fIJ.params->colidx[j];
        Bsym = fAB.params->qsym[b];
        Jsym = fIJ.params->qsym[j];
        if((I==J)&&(Asym==Bsym)) Amat.matrix[h][ai][bj] += fAB.matrix[Asym][A][B];
        if((A==B)&&(Isym==Jsym)) Amat.matrix[h][ai][bj] -= fIJ.matrix[Isym][I][J];
      }
    }
    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }
  dpd_buf4_close(&Amat);

  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fAB);
  dpd_file2_mat_close(&fIJ);
  dpd_file2_close(&fIJ);
}

void uhf_mohess(void)
{

}
