#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include "iwl.h"

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

/*
** iwl_rdtwo(): read two electron ints from the given file.
** The "iwl" stands for "integrals with labels," and this is the proposed
** new standard for storing two-electron integrals and their (absolute)
** orbital labels.
**
** Arguments:
**    itap     = unit to write to
**    ints     = two electron integrals (already allocated)
**    ioff     = the old ioff array for lexical ordering
**    norbs    = number of orbitals
**    nfzc     = number of frozen core orbitals
**    nfzv     = number of frozen virtual orbitals
**    printflg = print integrals as they're read 
**    outfile  = output file pointer
**
** David Sherrill, 1995
*/
void iwl_rdtwo(int itap, double *ints, int *ioff, int norbs, 
      int nfzc, int nfzv, int printflg, FILE *outfile)
{
  struct iwlbuf Buf;
  
  iwl_buf_init(&Buf, itap, 0.0, 1, 1);
  if ((nfzc == 0) && (nfzv == 0))
    iwl_buf_rd_all(&Buf, ints, ioff, ioff, 0, ioff, printflg, outfile);
  else
    iwl_buf_rd_all_act(&Buf, ints, ioff, ioff, 0, ioff, nfzc, norbs-nfzv-1,
                       printflg, outfile);
  iwl_buf_close(&Buf, 1);
}


