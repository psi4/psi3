#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include "iwl.h"

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))


/*
** iwl_buf_rd_all_mp2r12a()
**
** Read from an Integrals With Labels formatted buffer.
** The buffer must have been initialized with iwl_buf_init().
**
** Arguments:
**    Buf           =  IWL Buffer to read from (already initialized)
**    ints          =  memory buffer to put integrals into
**    ioff_lt       =  ioff array for the left pair of indices (p and q)
**    ioff_rt       =  ioff array for the right pair of indices (r and s)
**    bra_ket_symm  =  if 1, then these are ERI or R12 integrals, read
**                     them in as usual, else these are [r12,T2] integrals -
**                     form [T1+T2,r12] out of these.
**
**    WARNING - if bra_ket_symm = 0 - ints must be zeroed out!
**
**    ioff          =  the ioff array to figure the total index pqrs from
**                     the pair indices pq and rs
**    printflg      =  if 1, print integrals as they are read
**    outfile       =  pointer to output file for printing
**
** Returns: 0 if end of file, otherwise 1
**
*/
int iwl_buf_rd_all_mp2r12a(struct iwlbuf *Buf, double *ints,
			   int *ioff_lt, int *ioff_rt, int bra_ket_symm, int *ioff,
			   int printflg, FILE *outfile)
{
  int lastbuf;
  Label *lblptr;
  Value *valptr;
  int idx, p, q, r, s, pq, rs, pqrs;
  
  lblptr = Buf->labels;
  valptr = Buf->values;
  
  lastbuf = Buf->lastbuf;
  
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    pq = ioff_lt[p] + q;
    rs = ioff_rt[r] + s;

    pqrs = INDEX(pq,rs);

    if (bra_ket_symm) /* ERIs or R12-integrals */
      ints[pqrs] = (double) valptr[Buf->idx];
    else {            /* (ip|[T1+T2,r12]|jq) = - [ (ip|[r12,T1]|jq) + (jq|[r12,T2]|ip) ] */
      if (pq != rs)
	ints[pqrs] -= (double) valptr[Buf->idx];
      else
	ints[pqrs] -= (double) 2.0*valptr[Buf->idx];
    }
    
    if (printflg) 
      fprintf(outfile, "<%2d %2d %2d %2d [%2d][%2d] [[%3d]] = %20.10lf\n",
	      p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;
    
  } /* end loop through current buffer */
  
   /* read new PSI buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;
    
    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      pq = ioff_lt[p] + q;
      rs = ioff_rt[r] + s;

      pqrs = INDEX(pq,rs);

      if (bra_ket_symm) /* ERIs or R12-integrals */
	ints[pqrs] = (double) valptr[Buf->idx];
      else {            /* (ip|[T1+T2,r12]|jq) = - [ (ip|[r12,T2]|jq) + (jq|[r12,T2]|ip) ] */
	if (pq != rs)
	  ints[pqrs] -= (double) valptr[Buf->idx];
	else
	  ints[pqrs] -= (double) 2.0*valptr[Buf->idx];
      }
      
      if (printflg) 
	fprintf(outfile, "<%d %d %d %d [%d][%d] [[%d]] = %20.10lf\n",
		p, q, r, s, pq, rs, pqrs, ints[pqrs]) ;
      
    } /* end loop through current buffer */
    
  } /* end loop over reading buffers */
  
  return(0); /* we must have reached the last buffer at this point */
}

