#include <dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

/* fock_build(): Determine the contribution of a given integral to the
** molecular spin-orbital Fock matrices.  This routine assumes that
** the given integral indices are in canonical QT Standard ordering,
** i.e., p>=q, r>=s, and pq>=rs.  This convention is followed in
** transqt [see transform_two()] in the final call to libiwl's
** iwl_buf_wrt_mat().
**
** This algorithm is not currently used.  See fock.c.
**
** T. Daniel Crawford, February, 1999.  */

void fock_build(int p, int q, int r, int s, double value)
{
  int pq, rs, *dubble, *single, *virtual;
  double **fockA, **fockB;

  fockA = moinfo.fockA;
  fockB = moinfo.fockB;

  dubble = moinfo.dubble;
  single = moinfo.single;
  virtual = moinfo.virtual;
  
  pq = ioff[p] + q;  rs = ioff[r] + s;

  if(pq == rs) { /* p==r and q==s */
      if(r==s) { /* (pp|pp) */
	  if(dubble[r]) {
	      fockA[p][p] += value;
	      fockB[p][p] += value;
	    }
	  else if(single[r]) fockB[p][p] += value;
	}
      else { /* (pq|pq) p>q and r>s */
	  if(dubble[p] && dubble[q]) {
	      fockA[p][p] -= value;
	      fockB[p][p] -= value;
	      fockA[q][q] -= value;
	      fockB[q][q] -= value;
	    }
	  else if(single[p] && dubble[q]) {
	      fockA[q][q] -= value;
	      fockA[p][p] -= value;
	      fockB[p][p] -= value;
	    }
	  else if(single[p] && single[q]) {
	      fockA[p][p] -= value;
	      fockA[q][q] -= value;
	    }
	  else if(virtual[p] && dubble[q]) {
	      fockA[p][p] -= value;
	      fockB[p][p] -= value;
	    }
	  else if(virtual[p] && single[q]) {
	      fockA[p][p] -= value;
	    }
	}
    }
  else { /* p!=r &| q!=s */
      if(p==r) { /* q!=s */
	  if(p==q) { /* (pp|ps) */
	      if(dubble[p] && dubble[s]) {
		  fockA[p][s] += value;
		  fockB[p][s] += value;
		  fockA[s][p] += value;
		  fockB[s][p] += value;
		}
	      else if(single[p] && dubble[s]) {
		  fockB[p][s] += value;
		  fockB[s][p] += value;
		}
	      else if(single[p] && single[s]) {
		  fockB[p][s] += value;
		  fockB[s][p] += value;
		}
	    }
	  else { /* (pq|ps) */
	      if(dubble[p] && dubble[q] && dubble[s]) {
		  fockA[q][s] -= value;
		  fockA[s][q] -= value;
		  fockB[q][s] -= value;
		  fockB[s][q] -= value;
		}
	      else if(single[p] && dubble[q] && dubble[s]) {
		  fockA[q][s] -= value;
		  fockA[s][q] -= value;
		}
	      else if(single[p] && single[q] && dubble[s]) {
		  fockA[q][s] -= value;
		  fockA[s][q] -= value;
		}
	      else if(single[p] && single[q] && single[s]) {
		  fockA[q][s] -= value;
		  fockA[s][q] -= value;
		}
	    }
	}
      else if(q==s) { /* p!=r */
	  if(r==s) { /* (pq|qq) */
	      if(dubble[p] && dubble[q]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(single[p] && dubble[q]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(single[p] & single[q]) {
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(virtual[p] && dubble[q]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(virtual[p] && single[q]) {
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	    }
	  else { /* (pq|rq) */
	      if(dubble[p] && dubble[q] && dubble[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		  fockB[p][r] -= value;
		  fockB[r][p] -= value;
		}
	      else if(single[p] && dubble[q] && dubble[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		  fockB[p][r] -= value;
		  fockB[r][p] -= value;
		}
	      else if(single[p] && dubble[q] && single[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		  fockB[p][r] -= value;
		  fockB[r][p] -= value;
		}
	      else if(single[p] && single[q] && single[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		}
	      else if(virtual[p] && dubble[q] && dubble[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		  fockB[p][r] -= value;
		  fockB[r][p] -= value;
		}
	      else if(virtual[p] && dubble[q] && single[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		  fockB[p][r] -= value;
		  fockB[r][p] -= value;
		}
	      else if(virtual[p] && single[q] && single[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		}
	      else if(virtual[p] && dubble[q] && virtual[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		  fockB[p][r] -= value;
		  fockB[r][p] -= value;
		}
	      else if(virtual[p] && single[q] && virtual[r]) {
		  fockA[p][r] -= value;
		  fockA[r][p] -= value;
		}
	    }
	}
      else { /* p!=r && q!=s */
	  if(p==q && r==s) { /* (pp|rr) p>r */
	      if(dubble[p] && dubble[r]) {
		  fockA[p][p] += 2 * value;
		  fockB[p][p] += 2 * value;
		  fockA[r][r] += 2 * value;
		  fockB[r][r] += 2 * value;
		}
	      else if(single[p] && dubble[r]) {
		  fockA[r][r] += value;
		  fockB[r][r] += value;
		  fockA[p][p] += 2 * value;
		  fockB[p][p] += 2 * value;
		}
	      else if(single[p] && single[r]) {
		  fockA[p][p] += value;
		  fockA[r][r] += value;
		  fockB[p][p] += value;
		  fockB[r][r] += value;
		}
	      else if(virtual[p] && dubble[r]) {
		  fockA[p][p] += 2 * value;
		  fockB[p][p] += 2 * value;
		}
	      else if(virtual[p] && single[r]) {
		  fockA[p][p] += value;
		  fockB[p][p] += value;
		}
	    }
	  else if(p==q && r!=s) { /* (pp|rs) */
	      if(dubble[p] && dubble[r] && dubble[s]) {
		  fockA[r][s] += 2 * value;
		  fockA[s][r] += 2 * value;
		  fockB[r][s] += 2 * value;
		  fockB[s][r] += 2 * value;
		}
	      else if(single[p] && dubble[r] && dubble[s]) {
		  fockA[r][s] += value;
		  fockA[s][r] += value;
		  fockB[r][s] += value;
		  fockB[s][r] += value;
		}
	      else if(single[p] && single[r] && dubble[s]) {
		  fockA[r][s] += value;
		  fockA[s][r] += value;
		  fockB[r][s] += value;
		  fockB[s][r] += value;
		}
	      else if(single[p] && single[r] && single[s]) {
		  fockA[r][s] += value;
		  fockA[s][r] += value;
		  fockB[r][s] += value;
		  fockB[s][r] += value;
		}
	      
	    }
	  else if(p!=q && r==s) { /* (pq|rr) */
	      if(dubble[p] && dubble[q] && dubble[r]) {
		  fockA[p][q] += 2 * value;
		  fockA[q][p] += 2 * value;
		  fockB[p][q] += 2 * value;
		  fockB[q][p] += 2 * value;
		}
	      else if(single[p] && dubble[q] && dubble[r]) {
		  fockA[p][q] += 2 * value;
		  fockA[q][p] += 2 * value;
		  fockB[p][q] += 2 * value;
		  fockB[q][p] += 2 * value;
		}
	      else if(single[p] && dubble[q] && single[r]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(single[p] && single[q] && dubble[r]) {
		  fockA[p][q] += 2 * value;
		  fockA[q][p] += 2 * value;
		  fockB[p][q] += 2 * value;
		  fockB[q][p] += 2 * value;
		}
	      else if(single[p] && single[q] && single[r]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(virtual[p] && dubble[q] && dubble[r]) {
		  fockA[p][q] += 2 * value;
		  fockA[q][p] += 2 * value;
		  fockB[p][q] += 2 * value;
		  fockB[q][p] += 2 * value;
		}
	      else if(virtual[p] && dubble[q] && single[r]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(virtual[p] && single[q] && dubble[r]) {
		  fockA[p][q] += 2 * value;
		  fockA[q][p] += 2 * value;
		  fockB[p][q] += 2 * value;
		  fockB[q][p] += 2 * value;
		}
	      else if(virtual[p] && single[q] && single[r]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	      else if(virtual[p] && virtual[q] && dubble[r]) {
		  fockA[p][q] += 2 * value;
		  fockA[q][p] += 2 * value;
		  fockB[p][q] += 2 * value;
		  fockB[q][p] += 2 * value;
		}
	      else if(virtual[p] && virtual[q] && single[r]) {
		  fockA[p][q] += value;
		  fockA[q][p] += value;
		  fockB[p][q] += value;
		  fockB[q][p] += value;
		}
	    }
	  else if(p!=s && q==r) { /* (pq|qs) */
	      if(dubble[p] && dubble[q] && dubble[s]) {
		  fockA[p][s] -= value;
		  fockA[s][p] -= value;
		  fockB[p][s] -= value;
		  fockB[s][p] -= value;
		}
	      else if(single[p] && dubble[q] && dubble[s]) {
		  fockA[p][s] -= value;
		  fockA[s][p] -= value;
		  fockB[p][s] -= value;
		  fockB[s][p] -= value;
		}
	      else if(single[p] && single[q] && dubble[s]) {
		  fockA[p][s] -= value;
		  fockA[s][p] -= value;
		}
	      else if(single[p] && single[q] && single[s]) {
		  fockA[p][s] -= value;
		  fockA[s][p] -= value;
		}
	      else if(virtual[p] && dubble[q] && dubble[s]) {
		  fockA[p][s] -= value;
		  fockA[s][p] -= value;
		  fockB[p][s] -= value;
		  fockB[s][p] -= value;
		}
	      else if(virtual[p] && single[q] & dubble[s]) {
		  fockA[p][s] -= value;
		  fockA[s][p] -= value;
		}
	      else if(virtual[p] && single[q] && single[s]) {
		  fockA[p][s] -= value; 
		  fockA[s][p] -= value;
		}
	    }
	}
    }
}
