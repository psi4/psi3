/*****************************************************************************/
/* twoel_to_cas -  This routine writes out the two-electron integrals to a   */
/*                 form more palatable to the CASSCF code.  This should be   */
/*                 eliminated one day, and the CASSCF code should read the   */
/*                 two-electron integrals directly using libiwl.             */
/*                                                                           */
/*                 The two-electron integrals from TRANSQT are in CI order   */
/*                 but CASSCF wants Pitzer order.                            */
/*                 Use pqrs with p>=q, r>=s, pq>=rs                          */
/*                                                                           */
/* David Sherrill, April 1998                                                */
/*****************************************************************************/

#include <stdio.h>            /* the standard io library       */
#include <math.h>             /* the C math routines library   */
#include <libciomr.h>         /* standard psi routines library */
#include <qt.h>               /* the quantum trio library      */
#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)

/*
** declare needed external variables
*/
extern FILE *infile;
extern FILE *outfile;
extern int *ioff;

void twoel_to_cas(double *twoel_ints, int *corr_to_pitz, int nso,
                  int print_lvl, int itap)
{
  PSI_FPTR fptr = 0;            /* index for twoel ints for CAS output file */
  int ntri = 0;                 /* number of elements of lower triangle     */
  int cas_ints_dim = 0;         /* dimension of the twoel ints output buff  */
  double *cas_ints;             /* pointer to in-core twoel ints (CI order) */
  int p,q,r,s,pq,rs,smax;
  int so_p, so_q, so_r, so_s;
  int so_pq, so_rs;
  int pqrs,target;

  ntri = (nso * (nso+1))/2;
  cas_ints_dim = (ntri * (ntri+1))/2;

  cas_ints = init_array(cas_ints_dim);

  /* loop through all of the twoel ints and put in Pitzer order in buffer */
  for (p=0; p<nso; p++) {
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

          pq = ioff[p] + q;
          rs = ioff[r] + s;
          pqrs = INDEX(pq,rs);

          so_pq = INDEX(so_p, so_q);
          so_rs = INDEX(so_r, so_s);
          target = INDEX(so_pq, so_rs);

          cas_ints[target] = twoel_ints[pqrs];

        } 
      }
    }
  }


  /*
  ** write information to diskfile and close
  */
  rfile(itap);
  wwritw(itap, (char *) cas_ints, sizeof(double)*cas_ints_dim, 
         fptr, &fptr);  
  rclose(itap, 3);

  free(cas_ints);

}

