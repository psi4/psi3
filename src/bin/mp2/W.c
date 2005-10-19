#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void rhf_W(void);
void uhf_W(void);

void W(void)
{
  if(params.ref == 0) rhf_W();
  else if(params.ref == 2) uhf_W();
}

void rhf_W(void)
{
  int i, j, a, b, h;
  int nirreps;
  double fii = 0.0;
  double fjj = 0.0;
  double faa = 0.0;
  double fbb = 0.0;
  double Dij = 0.0;
  double Dab = 0.0;
  double Dai = 0.0;
  dpdfile2 W;
  dpdfile2 D;
  dpdfile2 f;
  dpdbuf4 T;
  dpdbuf4 I;

  nirreps = mo.nirreps;

  dpd_file2_init(&W, CC_OEI, 0, 0, 0, "WIJ");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&T, &I, &W, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_file2_close(&W);

  dpd_file2_init(&W, CC_OEI, 0, 1, 1, "WAB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&T, &I, &W, 3, 3, -2.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_file2_close(&W);

  dpd_file2_init(&W, CC_OEI, 0, 1, 0, "WAI");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_buf4_init(&I, CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
  dpd_contract442(&T, &I, &W, 2, 0, -4.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_file2_close(&W);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_file2_init(&f, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  dpd_file2_init(&W, CC_OEI, 0, 0, 0, "WIJ");
  dpd_file2_mat_init(&W);
  dpd_file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(i=0; i < mo.occpi[h]; i++) {
      fii = f.matrix[h][i][i];
      for(j=0; j < mo.occpi[h]; j++) {
        fjj = f.matrix[h][j][j];
        Dij = D.matrix[h][i][j];
        W.matrix[h][i][j] -= 0.5 * Dij * ( fii + fjj );
      }
    }
  }
  dpd_file2_mat_wrt(&W);
  dpd_file2_mat_close(&W);
  dpd_file2_close(&W);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_file2_init(&f, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  dpd_file2_init(&W, CC_OEI, 0, 1, 1, "WAB");
  dpd_file2_mat_init(&W);
  dpd_file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(a=0; a < mo.virpi[h]; a++) {
      faa = f.matrix[h][a][a];
      for(b=0; b < mo.virpi[h]; b++) {
        fbb = f.matrix[h][b][b];
        Dab = D.matrix[h][a][b];
        W.matrix[h][a][b] -= 0.5 * Dab * ( faa + fbb );
      }
    }
  }

  dpd_file2_mat_wrt(&W);
  dpd_file2_mat_close(&W);
  dpd_file2_close(&W);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_file2_init(&f, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  dpd_file2_init(&W, CC_OEI, 0, 1, 0, "WAI");
  dpd_file2_mat_init(&W);
  dpd_file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(a=0; a < mo.virpi[h]; a++) {
      for(i=0; i < mo.occpi[h]; i++) {
        fii = f.matrix[h][i][i];
        Dai = D.matrix[h][a][i];
        W.matrix[h][a][i] -= Dai * fii;
      }
    }
  }

  dpd_file2_mat_wrt(&W);
  dpd_file2_mat_close(&W);
  dpd_file2_close(&W);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);

  dpd_file2_init(&f, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  dpd_file2_init(&W, CC_OEI, 0, 0, 1, "WIA");
  dpd_file2_mat_init(&W);
  dpd_file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(a=0; a < mo.virpi[h]; a++) {
      for(i=0; i < mo.occpi[h]; i++) {
        faa = f.matrix[h][a][a];
        Dai = D.matrix[h][a][i];
        W.matrix[h][i][a] += Dai * faa;
      }
    }
  }

  dpd_file2_mat_wrt(&W);
  dpd_file2_mat_close(&W);
  dpd_file2_close(&W);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&W, CC_OEI, 0, 0, 0, "WIJ");
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&I, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_dot13(&D, &I, &W, 0, 0, -2.0, 1.0);
  dpd_dot14(&D, &I, &W, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&I, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_dot24(&D, &I, &W, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_init(&I, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
  dpd_dot23(&D, &I, &W, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&I);
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
  dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot13(&D, &I, &W, 0, 0, -2.0, 1.0);
  dpd_dot13(&D, &I, &W, 0, 1, -2.0, 1.0);
  dpd_dot14(&D, &I, &W, 0, 0, 1.0, 1.0);
  dpd_dot14(&D, &I, &W, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_file2_close(&D);
  dpd_file2_close(&W);

  dpd_file2_init(&f, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&f);
  dpd_file2_mat_rd(&f);

  dpd_file2_init(&W, CC_OEI, 0, 0, 0, "WIJ");
  dpd_file2_mat_init(&W);
  dpd_file2_mat_rd(&W);

  for(h=0; h < nirreps; h++) {
    for(i=0; i < mo.occpi[h]; i++) {
      fii = f.matrix[h][i][i];
      W.matrix[h][i][i] += 2.0 * fii;
    }
  }
  dpd_file2_mat_wrt(&W);
  dpd_file2_mat_close(&W);
  dpd_file2_close(&W);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);
}

void uhf_W(void)
{

}
