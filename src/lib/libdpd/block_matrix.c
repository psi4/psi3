/*
** dpd_block_matrix(): Allocates a contiguous block of memory for use
** as a 2-dimensional array.
**
** The memory is provided as a contiguous block of doubles, with an
** additional array of pointers used to mark the beginning of each
** row.  This allows transparent 2-dimensional-array style access, but
** keeps memory together such that it could be used in FORTRAN BLAS
** routines, for example.
**
** Prior to allocation, this routine checks the current status of
** dpd_default->memfree to make sure the malloc() request will not
** overrun the user-specified memory limits.  If there is insufficient
** memory available, entries are deleted from the dpd_file4_cache (in
** LRU order) until the memory limits are satisfied.  If, after
** deletion of the entire dpd_file4_cache (or at least until no other
** zero-priority entries remain), there is still insufficient memory
** available to satisfy the request, a NULL pointer is returned to the
** caller, indicating that either an out-of-core algorithm must be
** used, or the caller must exit().
**
** TDC, 6/24/00
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

double **dpd_block_matrix(int n, int m)
{
  int i;
  double **A, *B;

#ifdef DPD_TIMER
timer_on("block_mat");
#endif DPD_TIMER

  A = NULL;  B = NULL;

  while((dpd_default->memfree - n*m) < 0) {
      /* Delete cache entries until there's enough memory or no more cache */

      /* Priority-based cache */
      if(dpd_default->cachetype == 1) {
          if(dpd_file4_cache_del_low())
	      dpd_error("dpd_block_matrix: No memory left.", stderr);
	}

      /* Least-recently-used cache */
      else if(dpd_default->cachetype == 0) {
          if(dpd_file4_cache_del_lru())
	      dpd_error("dpd_block_matrix: No memory left.", stderr);
	}

      else dpd_error("LIBDPD Error: invalid cachetype.");
    }

  if(!m || !n) return(NULL);
  
  if((A = (double **) malloc(n * sizeof(double *)))==NULL) {
    fprintf(stderr,"dpd_block_matrix: trouble allocating memory \n");
    fprintf(stderr,"n = %d\n",n);
    exit(1);
    }

  if((B = (double *) malloc(m*n * sizeof(double)))==NULL) {
    fprintf(stderr,"dpd_block_matrix: trouble allocating memory \n");
    fprintf(stderr,"m = %d\n",m);
    exit(1);
    }

  bzero(B, m*n*sizeof(double));

  for (i = 0; i < n; i++) A[i] = &(B[i*m]);

  /* Decrement the global memory counter */
  dpd_default->memfree -= n*m;

#ifdef DPD_TIMER
timer_off("block_mat");
#endif DPD_TIMER

  return(A);
}

void dpd_free_block(double **array, int n, int m)
{
  if(array == NULL) return;
  free(array[0]);
  free(array);
  dpd_default->memfree += n*m;
}
