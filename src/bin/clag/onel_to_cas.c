/*****************************************************************************/
/* onel_to_cas -   This routine writes out the one-electron integrals to a   */
/*                 form more palatable to the CASSCF code.  This should be   */
/*                 eliminated one day, and the CASSCF code should read the   */
/*                 one-electron integrals directly using libiwl.             */
/*                                                                           */
/*                 The one-electron integrals from TRANSQT are in CI order   */
/*                 but CASSCF wants Pitzer order.                            */
/*                 Use pq with p>=q                                          */
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

void onel_to_cas(double *onel_ints, int *corr_to_pitz, int nso,
                  int print_lvl, int itap)
{
  PSI_FPTR fptr = 0;            /* index for onel ints for CAS output file  */
  int cas_ints_dim = 0;         /* dimension of the onel  ints output buff  */
  double *cas_ints;             /* pointer to in-core onel ints (CI order)  */
  int p,q,pq;
  int so_p, so_q, so_pq;

  cas_ints_dim = (nso * (nso+1))/2;
  cas_ints = init_array(cas_ints_dim);

  /* loop through all of the twoel ints and put in Pitzer order in buffer */
  for (p=0; p<nso; p++) {
    so_p = corr_to_pitz[p];

    for (q=0; q<=p; q++) {
      so_q = corr_to_pitz[q];

      pq = ioff[p] + q;
      so_pq = INDEX(so_p, so_q);

      cas_ints[so_pq] = onel_ints[pq];
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

