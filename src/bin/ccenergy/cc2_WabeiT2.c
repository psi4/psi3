#include <stdio.h>
#include <stdlib.h>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

void cc2_WabeiT2(void) {

  int rowx, colx, rowz, colz, ab;
  int GX, GZ, Ge, Gi, Gj, hxbuf, hzbuf;
  dpdfile2 t1, tia, tIA;
  dpdbuf4 Z, W, X;
  dpdbuf4 t2, t2a, t2b, tIJAB, tijab, tIjAb;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    /*     dpd_buf4_init(&X, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi"); */
    /*     dpd_buf4_scm(&X, 0); */
    /*     dpd_buf4_close(&X); */

    /** Begin out-of-core contract244 **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 0, 5, 0, 0, "ZAbIj");
    /*     dpd_buf4_init(&X, CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)"); */
    /*     dpd_buf4_sort_axpy(&X, CC2_HET1, rspq, 5, 11, "CC2 WAbEi", 1); */
    /*     dpd_buf4_close(&X); */
    dpd_buf4_init(&X, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    dpd_contract244(&t1, &X, &Z, 1, 2, 1, 1, 0);

    /* Symmetry Info */
    GX = X.file.my_irrep; 
    GZ = Z.file.my_irrep;

    dpd_file2_mat_init(&t1);
    dpd_file2_mat_rd(&t1);

    for(hxbuf=0; hxbuf < moinfo.nirreps; hxbuf++) {
      hzbuf = hxbuf;

      dpd_buf4_mat_irrep_row_init(&X, hxbuf);
      dpd_buf4_mat_irrep_row_init(&Z, hzbuf);

      /* Loop over rows of the X factor and the target */
      for(ab=0; ab < Z.params->rowtot[hzbuf]; ab++) {

	dpd_buf4_mat_irrep_row_zero(&X, hxbuf, ab);
	dpd_buf4_mat_irrep_row_rd(&X, hxbuf, ab);

/* 	dpd_buf4_mat_irrep_row_zero(&Z, hzbuf, ab); */

	for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	  Ge = Gi^hxbuf^GX;
	  Gj = Gi^hzbuf^GZ;

	  rowx = X.params->rpi[Ge];
	  colx = X.params->spi[Gi];
	  rowz = Z.params->rpi[Gi];
	  colz = Z.params->spi[Gj];

	  if(rowz && colz && rowx && colx) {
/* 	    	  C_DGEMM('n','n',rowz,colz,rowx,1.0, */
/* 	    		  &(t1.matrix[Gj][0][0]),rowx, */
/* 	    		  &(X.matrix[hxbuf][0][X.col_offset[hxbuf][Ge]]),colx,0.0, */
/* 	    		  &(Z.matrix[hzbuf][0][Z.col_offset[hzbuf][Gi]]),colz); */
	  }
	}

/* 	dpd_buf4_mat_irrep_row_wrt(&Z, hzbuf, ab); */
      }

      dpd_buf4_mat_irrep_row_close(&X, hxbuf);
      dpd_buf4_mat_irrep_row_close(&Z, hzbuf);
    }
    dpd_file2_mat_close(&t1);
    dpd_buf4_close(&X);
    /*   dpd_buf4_close(&Z); */

    dpd_buf4_sort_axpy(&Z, CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, srqp, 0, 5, "New tIjAb", 1);
    /** End out-of-core contract244 **/
    /*     dpd_buf4_close(&W); */

    /*     dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb"); */
    /*     dpd_buf4_axpy(&Z, &t2, 1); */
    /*     dpd_buf4_close(&t2); */
    /*     dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1); */
    dpd_buf4_close(&Z);

    dpd_file2_close(&t1);
  }

  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 7, 11, 7, 0, "CC2 WABEI (EI,A>B)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract244(&tIA, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2a, &tIJAB, 1);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);

    /*** BB ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 7, 11, 7, 0, "CC2 Wabei (ei,a>b)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract244(&tia, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    dpd_buf4_axpy(&t2a, &tijab, 1);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    dpd_contract244(&tIA, &W, &tIjAb, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZjIbA");
    dpd_buf4_init(&W, CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WaBeI (eI,aB)");
    dpd_contract244(&tia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

  }

  else if(params.ref == 2) { /*** UHF ***/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /*** AA ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 21, 7, 21, 7, 0, "CC2 WABEI (EI,A>B)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract244(&tIA, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_axpy(&t2a, &tIJAB, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIJAB);

    /*** BB ***/
    dpd_buf4_init(&W, CC2_HET1, 0, 31, 17, 31, 17, 0, "CC2 Wabei (ei,a>b)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_contract244(&tia, &W, &t2, 1, 0, 0, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 10, 17, "T (ji,a>b)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ji,a>b)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    dpd_buf4_axpy(&t2a, &tijab, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tijab);

    /*** AB ***/
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_init(&W, CC2_HET1, 0, 26, 28, 26, 28, 0, "CC2 WAbEi (Ei,Ab)");
    dpd_contract244(&tIA, &W, &tIjAb, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 29, 23, 29, 0, "CC2 ZjIbA");
    dpd_buf4_init(&W, CC2_HET1, 0, 25, 29, 25, 29, 0, "CC2 WaBeI (eI,aB)");
    dpd_contract244(&tia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 22, 28, "New tIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

  }

}
