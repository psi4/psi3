#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

static double 
d1diag_t1(int arrayid, char *arrayname)
{
  int h, nirreps, i;
  double **T, **C, *E, max;
  struct oe_dpdfile T1;

  nirreps = moinfo.nirreps;
  max = 0.0;

  dpd_oe_file_init(&T1, arrayid, 0, 1, arrayname, 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  for(h=0; h < nirreps; h++) {
      if(T1.params->rowtot[h] && T1.params->coltot[h]) {
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

static double
d1diag_subblock(double **Tave, int row0, int rown, int col0, int coln)
{
  int i,j;
  int nrow = rown - row0;
  int ncol = coln - col0;
  double max = 0.;
  double **Tsub;
  double **Tsq;
  double *E, **C;

  if (nrow && ncol) {
      Tsub  = block_matrix(nrow, ncol);
      Tsq  = block_matrix(nrow, nrow);

      for (i=row0; i<rown; i++) {
          for (j=col0; j<coln; j++) {
              Tsub[i-row0][j-col0] = Tave[i][j];
            }
        }

      newmm(Tsub, 0, Tsub, 1, Tsq, nrow, ncol, nrow, 1.0, 0.0);

      E = init_array(nrow);
      C = block_matrix(nrow, nrow);
      sq_rsp(nrow, nrow, Tsq, E, 0, C, 1e-12);

      /* Find maximum eigenvalue of T */
      for(i=0; i < nrow; i++) if(E[i] > max) max = E[i];

      free_block(C);
      free(E);
      free_block(Tsq);
      free_block(Tsub);
    }

  return max;
}

static double
d1diag_t1_rohf(int arrayid_a, char *arrayname_a,
               int arrayid_b, char *arrayname_b,
               int doprint)
{
  int h, nirreps, i, j;
  double **Tave, tmp, max;
  double max_ph=0.0, max_xp=0.0, max_hx=0.0;
  struct oe_dpdfile T1_a, T1_b;

  double max_ab, max_uhf;
  double t1diag=0.0;
  double max_all=0.0, max_anti=0.0, max_alpha=0.0, max_beta=0.0;
  int ncorr = 0;

  nirreps = moinfo.nirreps;

  dpd_oe_file_init(&T1_a, arrayid_a, 0, 1, arrayname_a, 0, outfile);
  dpd_oe_file_mat_init(&T1_a);
  dpd_oe_file_mat_rd(&T1_a, 0, outfile);

  dpd_oe_file_init(&T1_b, arrayid_b, 0, 1, arrayname_b, 0, outfile);
  dpd_oe_file_mat_init(&T1_b);
  dpd_oe_file_mat_rd(&T1_b, 0, outfile);

  for(h=0; h < nirreps; h++) {
      int nrow = T1_a.params->rowtot[h];
      int ncol = T1_a.params->coltot[h];
      int nopen = moinfo.openpi[h];
      ncorr += 2*(nrow-nopen) + nopen;
      if(nrow && ncol) {
         Tave = block_matrix(ncol, ncol);

         for (i=0; i<nrow; i++) {
             for (j=0; j<ncol; j++) {
                 Tave[i][j] = (T1_a.matrix[h][i][j] + T1_b.matrix[h][i][j])/2.;
                 if (i < nrow-nopen) {
                     if (j < ncol-nopen) {
                         t1diag += Tave[i][j]*Tave[i][j];
                       }
                     else {
                         t1diag += 2.0*Tave[i][j]*Tave[i][j];
                       }
                   }
                 else {
                     if (j < ncol-nopen) {
                         t1diag += 2.0*Tave[i][j]*Tave[i][j];
                       }
                   }
               }
           }

         tmp = d1diag_subblock(Tave, 0, nrow, 0, ncol);
         if (tmp > max_all) max_all = tmp;

         tmp = d1diag_subblock(Tave, 0, nrow-nopen, 0, ncol-nopen);
         if (tmp > max_ph) max_ph = tmp;

         tmp = d1diag_subblock(Tave, 0, nrow-nopen, ncol-nopen, ncol);
         if (tmp > max_hx) max_hx = tmp;

         tmp = d1diag_subblock(Tave, nrow-nopen, nrow, 0, ncol-nopen);
         if (tmp > max_xp) max_xp = tmp;
         for (i=0; i<nrow; i++) {
             for (j=0; j<ncol; j++) {
                 Tave[i][j] = (T1_a.matrix[h][i][j] - T1_b.matrix[h][i][j])/2.;
               }
           }

         tmp = d1diag_subblock(Tave, 0, nrow, 0, ncol);
         if (tmp > max_anti) max_anti = tmp;

         tmp = d1diag_subblock(T1_a.matrix[h], 0, nrow, 0, ncol);
         if (tmp > max_alpha) max_alpha = tmp;

         tmp = d1diag_subblock(T1_b.matrix[h], 0, nrow, 0, ncol);
         if (tmp > max_beta) max_beta = tmp;

         free_block(Tave);
       }
    }

  dpd_oe_file_mat_close(&T1_a);
  dpd_oe_file_close(&T1_a);

  dpd_oe_file_mat_close(&T1_b);
  dpd_oe_file_close(&T1_b);
  max_all = sqrt(max_all);
  max_ph = sqrt(max_ph);
  max_hx = sqrt(max_hx);
  max_xp = sqrt(max_xp);
  max_anti = sqrt(max_anti);
  max_alpha = sqrt(max_alpha);
  max_beta = sqrt(max_beta);
  t1diag = sqrt(t1diag/ncorr);

  max_uhf = (max_all>max_anti?max_all:max_anti);
  max_ab = (max_alpha>max_beta?max_alpha:max_beta);

  max = max_ph;
  if (max_hx > max) max = max_hx;
  if (max_xp > max) max = max_xp;

  if (doprint) {
      fprintf(outfile, "\tD(OCCSD) = %7.6f\n", max);
      fprintf(outfile, "\tD(ph)    = %7.6f\n", max_ph);
      fprintf(outfile, "\tD(hx)    = %7.6f\n", max_hx);
      fprintf(outfile, "\tD(xp)    = %7.6f\n", max_xp);
      fprintf(outfile, "\tD(ab)    = %7.6f\n", max_ab);
      fprintf(outfile, "\tD(a)     = %7.6f\n", max_alpha);
      fprintf(outfile, "\tD(b)     = %7.6f\n", max_beta);
      fprintf(outfile, "\tD(UCCSD) = %7.6f\n", max_uhf);
      fprintf(outfile, "\tD(E)     = %7.6f\n", max_all);
      fprintf(outfile, "\tD(E-)    = %7.6f\n", max_anti);
      fprintf(outfile, "\tT1Diag   = %7.6f\n", t1diag);

      fprintf(outfile, "\n");

      fprintf(outfile,
              "%5.4f&%5.4f&%5.4f&%5.4f&%5.4f&%5.4f&%5.4f&%5.4f&%5.4f&%5.4f&%5.4f\\\\\n",
              max, max_ph, max_hx, max_xp,
              max_ab, max_alpha, max_beta,
              max_uhf, max_all, max_anti,
              t1diag);
    }

  return max;
}

double d1diag(void)
{
  double norm = 0.0;
  if (moinfo.iopen) {
      norm = d1diag_t1_rohf(CC_OEI, "tia", CC_OEI, "tIA", 0);
    }
  else {
      norm = d1diag_t1(CC_OEI,"tIA");
    }
  return norm;
}

void d1diag_print(void)
{
  d1diag_t1_rohf(CC_OEI, "tia", CC_OEI, "tIA", 1);
}

