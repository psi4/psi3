/*! \file 
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* build_Z_ROHF():  Solve the orbital Z-vector equations for ROHF refs:
**
**    sum E,M A(AI,EM) D(orb)(E,M) = -X(A,I)
**
** where A(AI,EM) is the orbital Hessian computed in build_A(), X(A,I)
** is the orbital rotation gradient computed in build_X(), and
** D(orb)(E,M) is the final Z-vector we want. 
**
*/

void build_Z_ROHF(void)
{
  dpdbuf4 A;
  dpdfile2 X1, D;
  double **X, **T, **Y, **Z;
  int num_ai, h, nirreps, a, i, count, lastcol, rank, ai;
  int *virtpi, *occpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; openpi = moinfo.openpi; virtpi = moinfo.virtpi;

  /*** Construct the ai transformation matrix which places all singly
       occupied orbital combinations at the end of the vector ***/
  
  /* First compute the number of ai pairs */
  num_ai = 0;
  for(h=0; h < nirreps; h++)
    num_ai += virtpi[h] * occpi[h];

  /*  fprintf(outfile, "num_ai = %d\n", num_ai); */

  /* Malloc space for the transformation matrix */
  T = block_matrix(num_ai,num_ai);

  /* Now compute the row/column swaps we need and the number of zero
     columns*/
  for(h=0,count=0,rank=0; h < nirreps; h++)
    for(a=0; a < virtpi[h]; a++)
      for(i=0; i < occpi[h]; i++) {
	if((a >= (virtpi[h] - openpi[h])) &&
	   (i >= (occpi[h] - openpi[h])) )
	  T[count][count] = 0.0;
	else {
	  T[count][count] = 1.0;
	  rank++;
	}
	count++;
      }
  count = 0;
  lastcol = num_ai-1;
  for(h=0; h < nirreps; h++)
    for(a=0; a < virtpi[h]; a++)
      for(i=0; i < occpi[h] && lastcol > count; i++,count++) {
	if(T[count][count] == 0.0) {
	  while (T[lastcol][lastcol] == 0.0) lastcol--;
	  if(lastcol > count) {
	    T[count][lastcol] = T[lastcol][count] = 1.0;
	    T[lastcol][lastcol] = 0.0;
	  }
	}
      }

  /*  print_mat(T, num_ai, num_ai, outfile); */

  /*** Finished building the transformation matrix ***/

  /* Place all the elements of the orbital rotation gradient, X into a
     linear array, Z */
  dpd_file2_init(&X1, CC_MISC, 0, 1, 0, "X(A,I)");
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  num_ai = 0;
  for(h=0; h < nirreps; h++)
    num_ai += X1.params->rowtot[h]*X1.params->coltot[h];

  Z = block_matrix(1,num_ai);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < X1.params->rowtot[h]; a++)
      for(i=0; i < X1.params->coltot[h]; i++) 
	Z[0][count++] = -X1.matrix[h][a][i];


  dpd_file2_mat_close(&X1);
  dpd_file2_close(&X1);

  /* Push the zero elements of X to the end of the vector */
  X = block_matrix(1,num_ai);
  /*  newmm(T,0,Z,0,X,num_ai,num_ai,1,1.0,0.0); */
  if(num_ai)
    C_DGEMV('n',num_ai, num_ai, 1.0, &(T[0][0]), num_ai, &(Z[0][0]), 1, 0.0, &(X[0][0]), 1);

  /* Now, grab only irrep 0 of the orbital Hessian */
  dpd_buf4_init(&A, CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  dpd_buf4_mat_irrep_init(&A, 0);
  dpd_buf4_mat_irrep_rd(&A, 0);

  /* Move the zero rows and columns of the Hessian to the bottom.
     Note that as long as we won't be writing A back to disk, it's OK
     to put the product back in A.matrix[0]. */
  Y = block_matrix(num_ai,num_ai);
  /*  newmm(A.matrix[0],0,T,0,Y,num_ai,num_ai,num_ai,1.0,0.0); */
  if(num_ai) 
    C_DGEMM('n','n',num_ai,num_ai,num_ai,1.0,&(A.matrix[0][0][0]),num_ai,
	    &(T[0][0]),num_ai,0.0,&(Y[0][0]),num_ai);
  /*  newmm(T,0,Y,0,A.matrix[0],num_ai,num_ai,num_ai,1.0,0.0); */
  if(num_ai)
    C_DGEMM('n','n',num_ai,num_ai,num_ai,1.0,&(T[0][0]),num_ai,&(Y[0][0]),num_ai,
	    0.0,&(A.matrix[0][0][0]),num_ai);
  free_block(Y);

  /* Trying out Matt's Pople code --- way to go, Matt! */
  pople(A.matrix[0], X[0], rank, 1, 1e-12, outfile, 0);

  dpd_buf4_mat_irrep_close(&A, 0);
  dpd_buf4_close(&A);

  /* Now re-order the elements of X back to the DPD format */
  /*  newmm(T,0,X,0,Z,num_ai,num_ai,1,1.0,0.0); */
  if(num_ai) 
    C_DGEMV('n',num_ai, num_ai, 1.0, &(T[0][0]), num_ai, &(X[0][0]), 1, 0.0, &(Z[0][0]), 1);
  free_block(X);

  /* We don't need the transformation matrix anymore */
  free_block(T);

  /*
  for(ai=0; ai < num_ai; ai++) fprintf(outfile, "Z[%d] = %20.15f\n", ai, Z[0][ai]);
  */

  /* Build the orbital component of Dai --- we'll build these as separate
     spin cases for future simplicity (e.g., UHF-based codes)*/
  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++) {
	D.matrix[h][a][i] = Z[0][count++];

	if(a >= (virtpi[h] - openpi[h])) D.matrix[h][a][i] = 0.0;
      }
  dpd_file2_mat_wrt(&D);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++) 
      for(i=0; i < D.params->coltot[h]; i++) {
	D.matrix[h][a][i] = Z[0][count++];
	      
	if(i >= (occpi[h] - openpi[h])) D.matrix[h][a][i] = 0.0;
      }
  dpd_file2_mat_wrt(&D);
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* We're done with Z */
  free_block(Z);
}


