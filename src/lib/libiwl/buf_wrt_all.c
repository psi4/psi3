#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include "iwl.h"

/*
** iwl_buf_wrt_all()
**
** Write out two electron ints to IWL file.  Assume that the integrals
** are in ijkl canonical order (no spatial symmetry).
**
** Arguments:
**    itap     = unit to write to
**    nbfso    = number of basis functions in symmetry orbitals
**    ints     = two electron integrals 
**    ioff     = the old ioff array for lexical ordering
**    printflg = print flag (1 or 0)
**    outfile  =  output file
**
** David Sherrill, 6/27/96
**
*/
void iwl_buf_wrt_all(struct iwlbuf *Buf, int nbfso, double *ints, int *ioff,
      int printflg, FILE *outfile)
{
  int idx, i, p, q, r, s, smax, pq, rs, pqrs;
  Label *lblptr;
  Value *valptr;
  
  /* go through the lexical order and print to the output file */
  for (p=0; p<nbfso; p++) {
    for (q=0; q<=p; q++) {
      pq = ioff[p] + q;
      for (r=0; r<=p; r++) {
	smax = (p==r) ? (q+1) : (r+1);
	for (s=0; s < smax; s++) {
	  rs = ioff[r] + s;
	  pqrs = ioff[pq] + rs;
	  if (fabs(ints[pqrs]) > Buf->cutoff) {
	    idx = 4 * Buf->idx;
	    lblptr[idx++] = (Label) p;
	    lblptr[idx++] = (Label) q;
	    lblptr[idx++] = (Label) r;
	    lblptr[idx++] = (Label) s;
	    valptr[Buf->idx] = (Value) ints[pqrs];
	    Buf->idx++;
	    if (printflg) fprintf(outfile, "%d %d %d %d [%d] = %10.6lf\n",
				  p, q, r, s, pqrs, ints[pqrs]) ;
	    
	    if (Buf->idx == Buf->ints_per_buf) {
	      Buf->lastbuf = 0;
	      Buf->inbuf = Buf->idx;
	      iwl_buf_put(Buf);
	      Buf->idx = 0;
	    } 
	  }
	}
      }
    }
  }
}


