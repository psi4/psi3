#include <dpd.h>
#define EXTERN
#include "globals.h"

/* BUILD_A(): Construct the molecular orbital Hessian, A. At the
** moment were actually building all symmetry blocks of A, though for
** the orbital Z-vector equations we really only need the totally
** symmetric components.
** */

void build_A(void)
{
  int h, nirreps, e, m, a, i, em, ai, E, M, A, I;
  int Esym, Msym, Asym, Isym;
  int *virtpi, *openpi, *occpi, *occ_off, *vir_off;
  int *qt_occ, *qt_vir; /* Spatial orbital translators */
  struct oe_dpdfile fIJ, fij, fAB, fab, fIA, fia;
  struct dpdbuf Amat, Amat2, D, C;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; openpi = moinfo.openpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  qt_occ = moinfo.qt_occ; qt_vir = moinfo.qt_vir;

  /* Two-electron integral contributions */
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_sort(&D, CC_MISC, rpsq, 11, 11, "A(EM,AI)", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&Amat, CC_MISC, 11, 11, 11, 11, 0, "A(EM,AI)", 0, outfile);
  dpd_buf_sort(&Amat, CC_TMP0, psrq, 11, 11, "D <im|ea> (ei,am)", 0, outfile);
  dpd_scm(&Amat, 2.0, 0, outfile);
  dpd_copy(&Amat, CC_TMP0, "A(EM,ai)", 0, outfile);
  dpd_buf_init(&D, CC_TMP0, 11, 11, 11, 11, 0, "D <im|ea> (ei,am)", 0, outfile);
  dpd_axpy(&D, &Amat, -1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_sort(&C, CC_TMP0, qpsr, 11, 11, "C <ai|bj>", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&C, CC_TMP0, 11, 11, 11, 11, 0, "C <ai|bj>", 0, outfile);
  dpd_axpy(&C, &Amat, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_copy(&Amat, CC_TMP0, "A(em,ai)", 0, outfile);
  dpd_buf_close(&Amat);

  /* Fock matrix contributions */
  dpd_oe_file_init(&fIJ, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_file_mat_init(&fIJ);
  dpd_oe_file_mat_rd(&fIJ, 0, outfile);
  dpd_oe_file_init(&fij, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_file_mat_init(&fij);
  dpd_oe_file_mat_rd(&fij, 0, outfile);
  dpd_oe_file_init(&fAB, CC_OEI, 1, 1, "fAB", 0, outfile);
  dpd_oe_file_mat_init(&fAB);
  dpd_oe_file_mat_rd(&fAB, 0, outfile);
  dpd_oe_file_init(&fab, CC_OEI, 1, 1, "fab", 0, outfile);
  dpd_oe_file_mat_init(&fab);
  dpd_oe_file_mat_rd(&fab, 0, outfile);
  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_mat_init(&fIA);
  dpd_oe_file_mat_rd(&fIA, 0, outfile);
  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_mat_init(&fia);
  dpd_oe_file_mat_rd(&fia, 0, outfile);

  dpd_buf_init(&Amat, CC_MISC, 11, 11, 11, 11, 0, "A(EM,AI)", 0,
	       outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Amat, h);
      dpd_buf_mat_irrep_rd(&Amat, h, 0, outfile);

      for(em=0; em < Amat.params->rowtot[h]; em++) {
	  e = Amat.params->roworb[h][em][0];
	  m = Amat.params->roworb[h][em][1];
	  E = fAB.params->rowidx[e];
	  M = fIJ.params->rowidx[m];
	  Esym = fAB.params->psym[e];
	  Msym = fIJ.params->psym[m];
	  for(ai=0; ai < Amat.params->coltot[h]; ai++) {
	      a = Amat.params->colorb[h][ai][0];
	      i = Amat.params->colorb[h][ai][1];
	      A = fAB.params->colidx[a];
	      I = fIJ.params->colidx[i];
	      Asym = fAB.params->qsym[a];
	      Isym = fIJ.params->qsym[i];

	      if((M==I) && (Esym==Asym))
		  Amat.matrix[h][em][ai] += fAB.matrix[Esym][E][A];
	      if((E==A) && (Msym==Isym))
		  Amat.matrix[h][em][ai] -= fIJ.matrix[Msym][M][I];

	      /* Check to see if these virtual indices actually
		 correspond to open-shell orbitals --- if so, set this
		 element to zero */
	      if((E >= (virtpi[Esym] - openpi[Esym])) ||
		 (A >= (virtpi[Asym] - openpi[Asym])) )
		  Amat.matrix[h][em][ai] = 0.0;
	    }
	}

      dpd_buf_mat_irrep_wrt(&Amat, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Amat, h);
    }

  dpd_buf_close(&Amat);

  dpd_buf_init(&Amat, CC_TMP0, 11, 11, 11, 11, 0, "A(em,ai)", 0,
	       outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Amat, h);
      dpd_buf_mat_irrep_rd(&Amat, h, 0, outfile);

      for(em=0; em < Amat.params->rowtot[h]; em++) {
	  e = Amat.params->roworb[h][em][0];
	  m = Amat.params->roworb[h][em][1];
	  E = fab.params->rowidx[e];
	  M = fij.params->rowidx[m];
	  Esym = fab.params->psym[e];
	  Msym = fij.params->psym[m];
	  for(ai=0; ai < Amat.params->coltot[h]; ai++) {
	      a = Amat.params->colorb[h][ai][0];
	      i = Amat.params->colorb[h][ai][1];
	      A = fab.params->colidx[a];
	      I = fij.params->colidx[i];
	      Asym = fab.params->qsym[a];
	      Isym = fij.params->qsym[i];

	      if((M==I) && (Esym==Asym))
		  Amat.matrix[h][em][ai] += fab.matrix[Esym][E][A];
	      if((E==A) && (Msym==Isym))
		  Amat.matrix[h][em][ai] -= fij.matrix[Msym][M][I];

	      /* Check to see if these occupied indices actually
		 correspond to open-shell orbitals --- if so, set this
		 element to zero */
	      if((M >= (occpi[Msym] - openpi[Msym])) ||
		 (I >= (occpi[Isym] - openpi[Isym])) )
		  Amat.matrix[h][em][ai] = 0.0;
	    }
	}

      dpd_buf_mat_irrep_wrt(&Amat, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Amat, h);
    }

  dpd_buf_close(&Amat);

  dpd_buf_init(&Amat, CC_TMP0, 11, 11, 11, 11, 0, "A(EM,ai)", 0,
	       outfile);

  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Amat, h);
      dpd_buf_mat_irrep_rd(&Amat, h, 0, outfile);

      for(em=0; em < Amat.params->rowtot[h]; em++) {
	  e = Amat.params->roworb[h][em][0];
	  m = Amat.params->roworb[h][em][1];
	  Esym = Amat.params->psym[e];
	  Msym = Amat.params->qsym[m];
	  E = e - vir_off[Esym];
	  M = m - occ_off[Msym];
	  for(ai=0; ai < Amat.params->coltot[h]; ai++) {
	      a = Amat.params->colorb[h][ai][0];
	      i = Amat.params->colorb[h][ai][1];
	      Asym = Amat.params->rsym[a];
	      Isym = Amat.params->ssym[i];
	      A = a - vir_off[Asym];
	      I = i - occ_off[Isym];

	      /* This comparison is somewhat tricky.  The algebraic
		 expression for the Fock matrix shift here is:

		 A(EM,ai) += delta(M,a) f(E,i)(beta)

		 The Kronecker Delta is actually a comparison between
		 the *spatial* orbitals associated with M, and A.
		 Hence we have to compare the spatial orbital
		 translation of the the two absolute orbital indices. */
	      if((qt_occ[m] == qt_vir[a]) && (Esym==Isym))
		  Amat.matrix[h][em][ai] += fia.matrix[Isym][I][E];

	      /* Check to see if these occupied and virtual indices
		 actually correspond to open-shell orbitals --- if so,
		 set this element to zero */
	      if((E >= (virtpi[Esym] - openpi[Esym])) ||
		 (I >= (occpi[Isym] - openpi[Isym])) )
		  Amat.matrix[h][em][ai] = 0.0;
	    }
	}

      dpd_buf_mat_irrep_wrt(&Amat, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Amat, h);
    }
  dpd_buf_sort(&Amat, CC_TMP0, rspq, 11, 11, "A(em,AI)", 0, outfile);
  dpd_buf_close(&Amat);

  dpd_oe_file_mat_close(&fIJ);
  dpd_oe_file_close(&fIJ);
  dpd_oe_file_mat_close(&fij);
  dpd_oe_file_close(&fij);
  dpd_oe_file_mat_close(&fAB);
  dpd_oe_file_close(&fAB);
  dpd_oe_file_mat_close(&fab);
  dpd_oe_file_close(&fab);
  dpd_oe_file_mat_close(&fIA);
  dpd_oe_file_close(&fIA);
  dpd_oe_file_mat_close(&fia);
  dpd_oe_file_close(&fia);

  /* Now sum all three A-matrix components and divide by 2 */
  dpd_buf_init(&Amat, CC_MISC, 11, 11, 11, 11, 0, "A(EM,AI)", 0, outfile);
  dpd_buf_init(&Amat2, CC_TMP0, 11, 11, 11, 11, 0, "A(em,ai)", 0, outfile);
  dpd_axpy(&Amat2, &Amat, 1.0, 0, outfile);
  dpd_buf_close(&Amat2);
  dpd_buf_init(&Amat2, CC_TMP0, 11, 11, 11, 11, 0, "A(EM,ai)", 0, outfile);
  dpd_axpy(&Amat2, &Amat, 1.0, 0, outfile);
  dpd_buf_close(&Amat2);
  dpd_buf_init(&Amat2, CC_TMP0, 11, 11, 11, 11, 0, "A(em,AI)", 0, outfile);
  dpd_axpy(&Amat2, &Amat, 1.0, 0, outfile);
  dpd_buf_close(&Amat2);
  dpd_scm(&Amat, 0.5, 0, outfile);
  dpd_buf_close(&Amat);
}
