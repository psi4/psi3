#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void diis(int iter)
{
  int vector_file=90;  /* Error vector storage */
  int amp_file=91; /* Amplitude vector storage */
  int nvector=8;  /* Number of error vectors to keep */
  int h, nirreps;
  int row, col, word, p, q;
  int diis_cycle;
  int vector_length=0;
  struct oe_dpdfile T1, T1a, T1b;
  struct dpdbuf T2, T2a, T2b;
  PSI_FPTR junk;
  double *error;
  double **B, *C, **vector;
  double product, determinant, maximum;
  double start;

  nirreps = moinfo.nirreps;
  
  /* Compute the length of a single error vector */
  dpd_oe_file_init(&T1, CC_TMP0, 0, 1, "T1(I,A)", 0, outfile);
  dpd_buf_init(&T2a, CC_TMP0, 2, 7, 2, 7, 0, "T2(IJ,AB)", 0, outfile);
  dpd_buf_init(&T2b, CC_TMP0, 0, 5, 0, 5, 0, "T2(Ij,Ab)", 0, outfile);
  for(h=0; h < nirreps; h++) {
      vector_length += 2 * T1.params->rowtot[h] * T1.params->coltot[h];
      vector_length += 2 * T2a.params->rowtot[h] * T2a.params->coltot[h];
      vector_length += T2b.params->rowtot[h] * T2b.params->coltot[h];
    }
  dpd_oe_file_close(&T1);
  dpd_buf_close(&T2a);
  dpd_buf_close(&T2b);

  /* If we haven't already, open the vector files for reading/writing */
  if(iter == 1) { rfile(vector_file); rfile(amp_file); }

  /* Set the diis cycle value */
  diis_cycle = (iter-1) % nvector;

  /* Build the current error vector and dump it to disk */
  error = init_array(vector_length);
  word=0;
  dpd_oe_file_init(&T1a, CC_OEI, 0, 1, "New tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1a);
  dpd_oe_file_mat_rd(&T1a, 0, outfile);
  dpd_oe_file_init(&T1b, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1b);
  dpd_oe_file_mat_rd(&T1b, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	  for(col=0; col < T1a.params->coltot[h]; col++)
	      error[word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
  dpd_oe_file_mat_close(&T1a);
  dpd_oe_file_close(&T1a);
  dpd_oe_file_mat_close(&T1b);
  dpd_oe_file_close(&T1b);

  dpd_oe_file_init(&T1a, CC_OEI, 0, 1, "New tia", 0, outfile);
  dpd_oe_file_mat_init(&T1a);
  dpd_oe_file_mat_rd(&T1a, 0, outfile);
  dpd_oe_file_init(&T1b, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1b);
  dpd_oe_file_mat_rd(&T1b, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	  for(col=0; col < T1a.params->coltot[h]; col++)
	      error[word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
  dpd_oe_file_mat_close(&T1a);
  dpd_oe_file_close(&T1a);
  dpd_oe_file_mat_close(&T1b);
  dpd_oe_file_close(&T1b);
  
  dpd_buf_init(&T2a, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&T2b, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      dpd_buf_mat_irrep_rd(&T2a, h, 0, outfile);
      dpd_buf_mat_irrep_init(&T2b, h);
      dpd_buf_mat_irrep_rd(&T2b, h, 0, outfile);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf_mat_irrep_close(&T2a, h);
      dpd_buf_mat_irrep_close(&T2b, h);
    }
  dpd_buf_close(&T2a);
  dpd_buf_close(&T2b);

  dpd_buf_init(&T2a, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&T2b, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      dpd_buf_mat_irrep_rd(&T2a, h, 0, outfile);
      dpd_buf_mat_irrep_init(&T2b, h);
      dpd_buf_mat_irrep_rd(&T2b, h, 0, outfile);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf_mat_irrep_close(&T2a, h);
      dpd_buf_mat_irrep_close(&T2b, h);
    }
  dpd_buf_close(&T2a);
  dpd_buf_close(&T2b);

  dpd_buf_init(&T2a, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_buf_init(&T2b, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      dpd_buf_mat_irrep_rd(&T2a, h, 0, outfile);
      dpd_buf_mat_irrep_init(&T2b, h);
      dpd_buf_mat_irrep_rd(&T2b, h, 0, outfile);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf_mat_irrep_close(&T2a, h);
      dpd_buf_mat_irrep_close(&T2b, h);
    }
  dpd_buf_close(&T2a);
  dpd_buf_close(&T2b);

  start = diis_cycle*vector_length;
  wwritw(vector_file, (char *) error, vector_length*sizeof(double),
	 (unsigned long int) start*sizeof(double),&junk);

  /* Store the current amplitude vector on disk */
  word=0;
  dpd_oe_file_init(&T1a, CC_OEI, 0, 1, "New tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1a);
  dpd_oe_file_mat_rd(&T1a, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	  for(col=0; col < T1a.params->coltot[h]; col++)
	      error[word++] = T1a.matrix[h][row][col];
  dpd_oe_file_mat_close(&T1a);
  dpd_oe_file_close(&T1a);

  dpd_oe_file_init(&T1a, CC_OEI, 0, 1, "New tia", 0, outfile);
  dpd_oe_file_mat_init(&T1a);
  dpd_oe_file_mat_rd(&T1a, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	  for(col=0; col < T1a.params->coltot[h]; col++)
	      error[word++] = T1a.matrix[h][row][col];
  dpd_oe_file_mat_close(&T1a);
  dpd_oe_file_close(&T1a);
  
  dpd_buf_init(&T2a, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      dpd_buf_mat_irrep_rd(&T2a, h, 0, outfile);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      error[word++] = T2a.matrix[h][row][col];
      dpd_buf_mat_irrep_close(&T2a, h);
    }
  dpd_buf_close(&T2a);

  dpd_buf_init(&T2a, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      dpd_buf_mat_irrep_rd(&T2a, h, 0, outfile);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      error[word++] = T2a.matrix[h][row][col];
      dpd_buf_mat_irrep_close(&T2a, h);
    }
  dpd_buf_close(&T2a);

  dpd_buf_init(&T2a, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      dpd_buf_mat_irrep_rd(&T2a, h, 0, outfile);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      error[word++] = T2a.matrix[h][row][col];
      dpd_buf_mat_irrep_close(&T2a, h);
    }
  dpd_buf_close(&T2a);

  start = diis_cycle*vector_length;
  wwritw(amp_file, (char *) error, vector_length*sizeof(double),
	 (unsigned long int) start*sizeof(double),&junk);

  /* If we haven't run through enough iterations, set the correct dimensions
     for the extrapolation */
  if(!(iter >= (nvector))) {
      if(iter < 2) return; /* Leave if we can't extrapolate at all */
      nvector = iter;
    }

  /* Now grab the full set of error vectors from the file */
  vector = init_matrix(nvector, vector_length);
  for(p=0; p < nvector; p++) 
      wreadw(vector_file, (char *) vector[p], vector_length*sizeof(double),
	     (unsigned long int) p*vector_length*sizeof(double),&junk);

  /* Build B matrix of error vector products */
  B = init_matrix(nvector+1,nvector+1);

  for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++) {
	  dot_arr(vector[p], vector[q], vector_length, &product); 
	  B[p][q] = product;
	}

  for(p=0; p < nvector; p++) {
      B[p][nvector] = -1;
      B[nvector][p] = -1;
    }

  B[nvector][nvector] = 0;

  /* Find the maximum value in B and scale all its elements */
  maximum = fabs(B[0][0]);
  for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++)
	  if(fabs(B[p][q]) > maximum) maximum = fabs(B[p][q]);

  for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++)
	  B[p][q] /= maximum; 

  /* Build the constant vector */
  C = init_array(nvector+1);
  C[nvector] = -1;

  /* Solve the linear equations */
  flin(B, C, nvector+1, 1, &determinant);

  /* Grab the old amplitude vectors */
  for(p=0; p < nvector; p++)
      wreadw(amp_file, (char *) vector[p], vector_length*sizeof(double),
	     (unsigned long int) p*vector_length*sizeof(double), &junk);
  
  /* Build the new amplitude vector from the old ones */
  for(q=0; q < vector_length; q++) {
      error[q] = 0.0;
      for(p=0; p < nvector; p++)
	  error[q] += C[p] * vector[p][q];
    }

  /* Now place these elements into the DPD amplitude arrays */
  word=0;
  dpd_oe_file_init(&T1a, CC_OEI, 0, 1, "New tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1a);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	  for(col=0; col < T1a.params->coltot[h]; col++)
	      T1a.matrix[h][row][col] = error[word++];
  dpd_oe_file_mat_wrt(&T1a, 0, outfile);
  dpd_oe_file_mat_close(&T1a);
  dpd_oe_file_close(&T1a);

  dpd_oe_file_init(&T1a, CC_OEI, 0, 1, "New tia", 0, outfile);
  dpd_oe_file_mat_init(&T1a);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	  for(col=0; col < T1a.params->coltot[h]; col++)
	      T1a.matrix[h][row][col] = error[word++];
  dpd_oe_file_mat_wrt(&T1a, 0, outfile);
  dpd_oe_file_mat_close(&T1a);
  dpd_oe_file_close(&T1a);
  
  dpd_buf_init(&T2a, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      T2a.matrix[h][row][col] = error[word++];
      dpd_buf_mat_irrep_wrt(&T2a, h, 0, outfile);
      dpd_buf_mat_irrep_close(&T2a, h);
    }
  dpd_buf_close(&T2a);

  dpd_buf_init(&T2a, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      T2a.matrix[h][row][col] = error[word++];
      dpd_buf_mat_irrep_wrt(&T2a, h, 0, outfile);
      dpd_buf_mat_irrep_close(&T2a, h);
    }
  dpd_buf_close(&T2a);

  dpd_buf_init(&T2a, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	  for(col=0; col < T2a.params->coltot[h]; col++)
	      T2a.matrix[h][row][col] = error[word++];
      dpd_buf_mat_irrep_wrt(&T2a, h, 0, outfile);
      dpd_buf_mat_irrep_close(&T2a, h);
    }
  dpd_buf_close(&T2a);

  /* Release memory and return */
  free_matrix(vector, nvector);
  free_matrix(B, nvector+1);
  free(C);
  free(error);

  return;
}
