/*****************************************************************************/
/* twopdm_to_cas - This routine writes out the twopdm to a form more         */
/*                 palatable to the CASSCF code.  This should be eliminated  */
/*                 one day, and the CASSCF code should read the twopdm       */
/*                 directly.                                                 */
/*                                                                           */
/*                 The twopdm from DETCI has pqrs with pq>=rs, but we want   */
/*                 pqrs with p>=q, r>=s, pq>=rs for the CASSCF code.         */
/*                                                                           */
/*                 Only populated orbitals (frozen core or active orbitals,  */
/*                 i.e., everything but frozen virtuals) are relevant        */
/*                 but nevertheless write out full-dimension twopdm for now  */
/*                 Also, map to Pitzer order                                 */
/*                                                                           */
/* David Sherrill, March 1998                                                */
/*****************************************************************************/

#include <stdio.h>            /* the standard io library       */
#include <math.h>             /* the C math routines library   */
#include <libciomr/libciomr.h>         /* standard psi routines library */
#include <qt.h>               /* the quantum trio library      */
#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)

/*
** declare needed external variables
*/
extern FILE *infile;
extern FILE *outfile;
extern int *ioff;

void twopdm_to_cas(double *tpdm, int *corr_to_pitz, int nso,
                   int npop, int print_lvl, int itap)
{
  PSI_FPTR tpdm_fptr = 0;       /* index for two-pdm for CAS output file */
  int ntri = 0;                 /* number of elements of lower triangle  */
  int symm_tpdm_dim = 0;        /* dimension of the symmetrized twopdm   */
  double *symm_tpdm;            /* pointer to in core symmetrized twopdm */
  int p,q,r,s,pq,qp,rs,sr,smax;
  int so_p, so_q, so_r, so_s;
  int so_pq, so_rs;
  int pqrs,qprs,pqsr,qpsr,target;
  double value;

  ntri = (nso * (nso+1))/2;
  symm_tpdm_dim = (ntri * (ntri+1))/2;

  symm_tpdm = init_array(symm_tpdm_dim);

  /* loop through all of the twopdm and put it in the symmetrized version */
  for (p=0; p<npop; p++) {
    so_p = corr_to_pitz[p];

    for (q=0; q<=p; q++) {
      so_q = corr_to_pitz[q];

      for (r=0; r<=p; r++) {
        so_r = corr_to_pitz[r];

        if (r==p) smax = q+1;
        else smax = r+1;

        /* adding orbsym here would be smart...save some time...but this
         * is just a hack code anyway right now, so... */
        
        /* pqrsym = orbsym[p] ^ orbsym[q] ^ orbsym[r]; */

        for (s=0; s<smax; s++) {
          so_s = corr_to_pitz[s];

          /* if ((orbsym[s] ^ pqrsym) != 0) continue; */

          pq = p * npop + q;
          qp = q * npop + p;
          rs = r * npop + s;
          sr = s * npop + r;
          pqrs = INDEX(pq,rs);
          qprs = INDEX(qp,rs);
          pqsr = INDEX(pq,sr);
          qpsr = INDEX(qp,sr);

          so_pq = INDEX(so_p, so_q);
          so_rs = INDEX(so_r, so_s);
          target = INDEX(so_pq, so_rs);

          value = 0.25 * (tpdm[pqrs] + tpdm[qprs] + 
                          tpdm[pqsr] + tpdm[qpsr]);
          symm_tpdm[target] = value;

        } 
      }
    }
  }


  /*
  ** write information to diskfile and close
  */
  rfile(itap);
  wwritw(itap, (char *) symm_tpdm, sizeof(double)*symm_tpdm_dim, 
         tpdm_fptr, &tpdm_fptr);  
  rclose(itap, 3);

  free(symm_tpdm);

}

