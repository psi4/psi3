#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr.h>
#include <psio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void diis(int iter)
{
  int nvector=8;  /* Number of error vectors to keep */
  int h, nirreps;
  int row, col, word, p, q;
  int diis_cycle;
  int vector_length=0;
  dpdfile2 T1, T1a, T1b;
  dpdbuf4 T2, T2a, T2b;
  psio_address start, end, next;
  double *error;
  double **B, *C, **vector;
  double product, determinant, maximum;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    /* Compute the length of a single error vector */
    dpd_file2_init(&T1, CC_TAMPS, 0, 0, 1, "tIA");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      vector_length += T1.params->rowtot[h] * T1.params->coltot[h];
      vector_length += T2.params->rowtot[h] * T2.params->coltot[h];
    }
    dpd_file2_close(&T1);
    dpd_buf4_close(&T2);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = init_array(vector_length);

    word=0;
    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    dpd_file2_init(&T1b, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1b);
    dpd_file2_mat_rd(&T1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  error[word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);
    dpd_file2_mat_close(&T1b);
    dpd_file2_close(&T1b);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&T2b, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      dpd_buf4_mat_irrep_init(&T2b, h);
      dpd_buf4_mat_irrep_rd(&T2b, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2b, h);
    }
    dpd_buf4_close(&T2a);
    dpd_buf4_close(&T2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_ERR, "DIIS Error Vectors" , (char *) error, 
	       vector_length*sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word=0;

    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  error[word++] = T1a.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_AMP, "DIIS Amplitude Vectors" , (char *) error, 
	       vector_length*sizeof(double), start, &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if(!(iter >= (nvector))) {
      if(iter < 2) { /* Leave if we can't extrapolate at all */
	free(error);
	return; 
      }
      nvector = iter;
    }

    /* Now grab the full set of error vectors from the file */
    vector = init_matrix(nvector, vector_length);
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++) 
      psio_read(CC_DIIS_ERR, "DIIS Error Vectors", (char *) vector[p], 
		vector_length*sizeof(double), next, &next);

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
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++) 
      psio_read(CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *) vector[p], 
		vector_length*sizeof(double), next, &next);
  
    /* Build the new amplitude vector from the old ones */
    for(q=0; q < vector_length; q++) {
      error[q] = 0.0;
      for(p=0; p < nvector; p++)
	error[q] += C[p] * vector[p][q];
    }

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  T1a.matrix[h][row][col] = error[word++];
    dpd_file2_mat_wrt(&T1a);
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  T2a.matrix[h][row][col] = error[word++];
      dpd_buf4_mat_irrep_wrt(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    /* Release memory and return */
    free_matrix(vector, nvector);
    free_matrix(B, nvector+1);
    free(C);
    free(error);
  }
  else { /** ROHF, UHF **/
  
    /* Compute the length of a single error vector */
    dpd_file2_init(&T1, CC_TMP0, 0, 0, 1, "tIA");
    dpd_buf4_init(&T2a, CC_TMP0, 0, 2, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&T2b, CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      vector_length += 2 * T1.params->rowtot[h] * T1.params->coltot[h];
      vector_length += 2 * T2a.params->rowtot[h] * T2a.params->coltot[h];
      vector_length += T2b.params->rowtot[h] * T2b.params->coltot[h];
    }
    dpd_file2_close(&T1);
    dpd_buf4_close(&T2a);
    dpd_buf4_close(&T2b);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = init_array(vector_length);
    word=0;
    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    dpd_file2_init(&T1b, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1b);
    dpd_file2_mat_rd(&T1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  error[word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);
    dpd_file2_mat_close(&T1b);
    dpd_file2_close(&T1b);

    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    dpd_file2_init(&T1b, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_mat_init(&T1b);
    dpd_file2_mat_rd(&T1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  error[word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);
    dpd_file2_mat_close(&T1b);
    dpd_file2_close(&T1b);
  
    dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&T2b, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      dpd_buf4_mat_irrep_init(&T2b, h);
      dpd_buf4_mat_irrep_rd(&T2b, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2b, h);
    }
    dpd_buf4_close(&T2a);
    dpd_buf4_close(&T2b);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&T2b, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      dpd_buf4_mat_irrep_init(&T2b, h);
      dpd_buf4_mat_irrep_rd(&T2b, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2b, h);
    }
    dpd_buf4_close(&T2a);
    dpd_buf4_close(&T2b);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&T2b, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      dpd_buf4_mat_irrep_init(&T2b, h);
      dpd_buf4_mat_irrep_rd(&T2b, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2b, h);
    }
    dpd_buf4_close(&T2a);
    dpd_buf4_close(&T2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_ERR, "DIIS Error Vectors" , (char *) error, 
	       vector_length*sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word=0;
    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  error[word++] = T1a.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);

    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  error[word++] = T1a.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);
  
    dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  error[word++] = T2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_AMP, "DIIS Amplitude Vectors" , (char *) error, 
	       vector_length*sizeof(double), start, &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if(!(iter >= (nvector))) {
      if(iter < 2) { /* Leave if we can't extrapolate at all */
	free(error);
	return; 
      }
      nvector = iter;
    }

    /* Now grab the full set of error vectors from the file */
    vector = init_matrix(nvector, vector_length);
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++) 
      psio_read(CC_DIIS_ERR, "DIIS Error Vectors", (char *) vector[p], 
		vector_length*sizeof(double), next, &next);

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
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++) 
      psio_read(CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *) vector[p], 
		vector_length*sizeof(double), next, &next);
  
    /* Build the new amplitude vector from the old ones */
    for(q=0; q < vector_length; q++) {
      error[q] = 0.0;
      for(p=0; p < nvector; p++)
	error[q] += C[p] * vector[p][q];
    }

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_mat_init(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  T1a.matrix[h][row][col] = error[word++];
    dpd_file2_mat_wrt(&T1a);
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);

    dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_mat_init(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
	for(col=0; col < T1a.params->coltot[h]; col++)
	  T1a.matrix[h][row][col] = error[word++];
    dpd_file2_mat_wrt(&T1a);
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);
  
    dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  T2a.matrix[h][row][col] = error[word++];
      dpd_buf4_mat_irrep_wrt(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  T2a.matrix[h][row][col] = error[word++];
      dpd_buf4_mat_irrep_wrt(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    dpd_buf4_init(&T2a, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
	for(col=0; col < T2a.params->coltot[h]; col++)
	  T2a.matrix[h][row][col] = error[word++];
      dpd_buf4_mat_irrep_wrt(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    /* Release memory and return */
    free_matrix(vector, nvector);
    free_matrix(B, nvector+1);
    free(C);
    free(error);
  }

  return;
}
