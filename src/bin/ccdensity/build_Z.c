#include <dpd.h>
#include <qt.h>
#include <libciomr.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* build_Z():  Solve the orbital Z-vector equations:
**
**    sum E,M A(AI,EM) D(orb)(E,M) = -X(A,I)
**
** where A(AI,EM) is the orbital Hessian computed in build_A(), X(A,I)
** is the orbital rotation gradient computed in build_X(), and D(orb)(E,M)
** is the final Z-vector we want.
** */

void build_Z(void)
{
  struct dpdbuf A;
  struct oe_dpdfile X1, D;
  double **X, **T, **Y, **Z;
  int num_ai, h, nirreps, a, i, count, lastcol, rank;
  int *virtpi, *occpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; openpi = moinfo.openpi; virtpi = moinfo.virtpi;

  /*** Construct the ai transformation matrix which places all singly
    occupied orbital combinations at the end of the vector ***/
  
  /* First compute the number of ai pairs */
  num_ai = 0;
  for(h=0; h < nirreps; h++)
      num_ai += virtpi[h] * occpi[h];

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

  /*** Finished building the transformation matrix ***/

  /* Place all the elements of the orbital rotation gradient, X into a
     linear array, Z */
  dpd_oe_file_init(&X1, CC_MISC, 1, 0, "X(A,I)", 0, outfile);
  dpd_oe_file_mat_init(&X1);
  dpd_oe_file_mat_rd(&X1, 0, outfile);
  num_ai = 0;
  for(h=0; h < nirreps; h++)
      num_ai += X1.params->rowtot[h]*X1.params->coltot[h];

  Z = block_matrix(1,num_ai);
  for(h=0,count=0; h < nirreps; h++)
      for(a=0; a < X1.params->rowtot[h]; a++)
	  for(i=0; i < X1.params->coltot[h]; i++)
	      Z[0][count++] = -X1.matrix[h][a][i];

  dpd_oe_file_mat_close(&X1);
  dpd_oe_file_close(&X1);

  /* Push the zero elements of X to the end of the vector */
  X = block_matrix(1,num_ai);
  newmm(T,0,Z,0,X,num_ai,num_ai,1,1.0,0.0);

  /* Now, grab only irrep 0 of the orbital Hessian */
  dpd_buf_init(&A, CC_MISC, 11, 11, 11, 11, 0, "A(EM,AI)", 0, outfile);
  dpd_buf_mat_irrep_init(&A, 0);
  dpd_buf_mat_irrep_rd(&A, 0, 0, outfile);

  /* Move the zero rows and columns of the Hessian to the bottom.
     Note that as long as we won't be writing A back to disk, it's OK
     to put the product back in A.matrix[0]. */
  Y = block_matrix(num_ai,num_ai);
  newmm(A.matrix[0],0,T,0,Y,num_ai,num_ai,num_ai,1.0,0.0);
  newmm(T,0,Y,0,A.matrix[0],num_ai,num_ai,num_ai,1.0,0.0);
  free_block(Y);

  /* Trying out Matt's Pople code --- way to go, Matt! */
  pople(A.matrix[0], X[0], rank, 1, 1e-12, outfile, 0);

  dpd_buf_mat_irrep_close(&A, 0);
  dpd_buf_close(&A);

  /* Now re-order the elements of X back to the DPD format */
  newmm(T,0,X,0,Z,num_ai,num_ai,1,1.0,0.0);
  free_block(X);

  /* We don't need the transformation matrix anymore */
  free_block(T);

  /* Build the orbital component of Dai --- we'll build these as separate
     spin cases for future simplicity (e.g., UHF-based codes)*/
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(A,I)", 0, outfile);
  dpd_oe_file_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
      for(a=0; a < D.params->rowtot[h]; a++)
	  for(i=0; i < D.params->coltot[h]; i++) {
	      D.matrix[h][a][i] = Z[0][count++];

	      if(a >= (virtpi[h] - openpi[h])) D.matrix[h][a][i] = 0.0;
	    }
  dpd_oe_file_mat_wrt(&D, 0, outfile);
  dpd_oe_file_mat_close(&D);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(a,i)", 0, outfile);
  dpd_oe_file_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
      for(a=0; a < D.params->rowtot[h]; a++) 
	  for(i=0; i < D.params->coltot[h]; i++) {
	      D.matrix[h][a][i] = Z[0][count++];
	      
	      if(i >= (occpi[h] - openpi[h])) D.matrix[h][a][i] = 0.0;
	    }
  dpd_oe_file_mat_wrt(&D, 0, outfile);
  dpd_oe_file_mat_close(&D);
  dpd_oe_file_close(&D);

  /* We're done with Z */
  free_block(Z);
}


