/*!
  \file
*/
#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include "iwl.h"
#include "iwl.hpp"

using namespace psi;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
      
int IWL::read_array2(double *ints, int *plist, int *qlist, int *rlist, int *slist,
    int *size, int *ioff, int printflg, FILE *outfile)
{
    int lastbuf;
    Label *lblptr;
    Value *valptr;
    int idx, p, q, pq, r, s;
    double value;

    lblptr = labels_;
    valptr = values_;

    lastbuf = lastbuf_;

    *size = 0;

    for (idx=4*idx_; idx_ < inbuf_; idx_++) {
        p = (int) lblptr[idx++];
        q = (int) lblptr[idx++];
        r = (int) lblptr[idx++];
        s = (int) lblptr[idx++];

        pq = INDEX(p,q);

        value = (double) valptr[idx_];
        *ints++ = value;
        *plist++ = p;
        *qlist++ = q;
        *rlist++ = r;
        *slist++ = s;
        *size= *size + 1;

        if (printflg) 
            fprintf(outfile, "<%d %d %d %d [%d] = %20.10f\n", p, q, r, s,
            pq, value);
    } /*! end loop through current buffer */

    /*! read new buffers */
    while (!lastbuf) {
        fetch();
        lastbuf = lastbuf_;

        for (idx=4*idx_; idx_ < inbuf_; idx_++) {
            p = (int) lblptr[idx++];
            q = (int) lblptr[idx++];
            r = (int) lblptr[idx++];
            s = (int) lblptr[idx++];

            pq = INDEX(p,q);

            value = (double) valptr[idx_];
            *ints++ = value;
            *plist++ = p;
            *qlist++ = q;
            *rlist++ = r;
            *slist++ = s;
            *size = *size + 1;

            if (printflg) 
                fprintf(outfile, "<%d %d %d %d [%d] = %20.10f\n", p, q, r, s,
                pq, value);
        } /*! end loop through current buffer */
    } /*! end loop over reading buffers */

    return(0); /*! we must have reached the last buffer at this point */
}

extern "C" {

/*!
** iwl_buf_rd_arr2()
**
** Read from an Integrals With Labels formatted PSI buffer.
** The buffer must have been initialized with iwl_buf_init().  The
** integrals and their labels are returned in the arrays ints, plist,
** qlist, rlist, and slist, and the size of these arrays is returned in 
** 'size.'
**
** Returns: 0 if end of file, otherwise 1
**
** Revised 6/27/96 by CDS for new format
**
*/
int iwl_buf_rd_arr2(struct iwlbuf *Buf, double *ints, int *plist, 
      int *qlist, int *rlist, int *slist, int *size, int *ioff,
      int printflg, FILE *outfile)
{
  int lastbuf;
  Label *lblptr;
  Value *valptr;
  int idx, p, q, pq, r, s;
  double value;
  
  lblptr = Buf->labels;
  valptr = Buf->values;
  
  lastbuf = Buf->lastbuf;
  
  *size = 0;
  
  for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];
    
    pq = INDEX(p,q);
    
    value = (double) valptr[Buf->idx];
    *ints++ = value;
    *plist++ = p;
    *qlist++ = q;
    *rlist++ = r;
    *slist++ = s;
    *size= *size + 1;
    
    if (printflg) 
      fprintf(outfile, "<%d %d %d %d [%d] = %20.10lf\n", p, q, r, s,
	      pq, value);
    
  } /*! end loop through current buffer */

  /*! read new buffers */
  while (!lastbuf) {
    iwl_buf_fetch(Buf);
    lastbuf = Buf->lastbuf;

    for (idx=4*Buf->idx; Buf->idx<Buf->inbuf; Buf->idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];
      
      pq = INDEX(p,q);

      value = (double) valptr[Buf->idx];
      *ints++ = value;
      *plist++ = p;
      *qlist++ = q;
      *rlist++ = r;
      *slist++ = s;
      *size = *size + 1;
      
      if (printflg) 
	fprintf(outfile, "<%d %d %d %d [%d] = %20.10lf\n", p, q, r, s,
		pq, value);
      
    } /*! end loop through current buffer */
    
  } /*! end loop over reading buffers */

  return(0); /*! we must have reached the last buffer at this point */
}

} /* extern "C" */