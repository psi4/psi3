#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))

void diag_A(void)
{
  int nirreps, h, h0, h1, a, i, j, *num_ai, count, lastcol, *rank, errcod;
  int *virtpi, *occpi, *openpi, ***dp;
  char *ref;
  double **T, **Y, *X, **evals;
  dpdbuf4 A;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  openpi = moinfo.openpi;
  virtpi = moinfo.virtpi;

  /* Final list of evals for pretty output */
  evals = init_matrix(nirreps,5);

  dpd_buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  /* Construct binary direct product array */
  dp = (int ***) malloc(nirreps * sizeof(int **));
  num_ai = init_int_array(nirreps);
  rank = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
      dp[h] = init_int_matrix(nirreps,2);
      count=0;
      for(h0=0; h0 < nirreps; h0++) {
          for(h1=0; h1 < nirreps; h1++) {
              if((h0^h1)==h) {
                  dp[h][count][0] = h0;
                  dp[h][count++][1] = h1;

		  num_ai[h] += virtpi[h0] * occpi[h1];
                }
            }
        }
    }


  for(h=0; h < nirreps; h++) {

      /*** Construct the ai transformation matrix which places all singly
	occupied orbital combinations at the end of the vector ***/
  
      /* Malloc space for the transformation matrix */
      T = block_matrix(num_ai[h],num_ai[h]);

      /* Now compute the row/column swaps we need and the number
	 of zero columns*/
      count = 0;
      rank[h] = 0;
      for(j=0; j < nirreps; j++) {
	  h0 = dp[h][j][0];  h1 = dp[h][j][1];
	  
	  for(a=0; a < virtpi[h0]; a++)
	      for(i=0; i < occpi[h1]; i++) {
		  
		  if((a >= (virtpi[h0] - openpi[h0])) &&
		     (i >= (occpi[h1] - openpi[h1])) )
		      T[count][count] = 0.0;
		  else {
		      T[count][count] = 1.0;
		      rank[h]++;
		    }
		  count++;
		}
	}
      
      count = 0;
      lastcol = num_ai[h]-1;
      for(j=0; j < nirreps; j++) {
	  h0 = dp[h][j][0];  h1 = dp[h][j][1];
	  
	  for(a=0; a < virtpi[h0]; a++)
	      for(i=0; i < occpi[h1] && lastcol > count; i++,count++) {
		  if(T[count][count] == 0.0) {
		      while (T[lastcol][lastcol] == 0.0) lastcol--;
		      if(lastcol > count) {
			  T[count][lastcol] = T[lastcol][count] = 1.0;
			  T[lastcol][lastcol] = 0.0;
			}
		    }
		}
	}

      /*** Finished building the transformation matrix ***/

      /* Apply T to move the zero rows and cols of the Hessian to the
	 bottom */
      dpd_buf4_mat_irrep_init(&A, h);
      dpd_buf4_mat_irrep_rd(&A, h);

      Y = block_matrix(num_ai[h], num_ai[h]); /* Scratch array */
      newmm(A.matrix[h], 0, T, 0, Y, num_ai[h], num_ai[h], num_ai[h], 1.0, 0.0);
      newmm(T, 0, Y, 0, A.matrix[h], num_ai[h], num_ai[h], num_ai[h], 1.0, 0.0);

      /* Get rid of this T */
      free_block(T);

      /* Diagonalize A --- does this need to be OOC?? */
      X = init_array(rank[h]);
      sq_rsp(rank[h], rank[h], A.matrix[h], X, 1, Y, 1e-14);

      for(i=0; i < MIN0(rank[h],5); i++)
	  evals[h][i] = X[i];

      free(X);
      free_block(Y);
    }

  /* Print the eigenvalues in a nice format */
  errcod = ip_string("REFERENCE",&ref,0);
  fprintf(outfile, "\n    Lowest eigenvalues of the %s orbital Hessian:\n",
	  ref);
  fprintf(outfile, "\n");
  free(ref);

  for(h=0; h < nirreps; h++)
      fprintf(outfile, "    %3s  ",moinfo.labels[h]);
  fprintf(outfile, "\n");
  for(h=0; h < nirreps; h++)
      fprintf(outfile, "---------");
  fprintf(outfile, "\n");

  for(i=0; i < 5; i++) {
      for(h=0; h < nirreps; h++) {
	  if(rank[h] <= i) fprintf(outfile, "         ");
	  else fprintf(outfile, " %7.4f ", evals[h][i]);
	}
      fprintf(outfile, "\n");
    }

  fprintf(outfile, "\n");

  free_matrix(evals, nirreps);
  free(rank);

  /* Get rid of num_ai and dp */
  free(num_ai);

  for(h=0; h < nirreps; h++)
      free_int_matrix(dp[h], nirreps);
  free(dp);

  dpd_buf4_close(&A);
}

