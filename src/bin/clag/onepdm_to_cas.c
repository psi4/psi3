/*****************************************************************************/
/* opdm_to_cas -   This routine writes out the one-particle density matrix   */
/*                 to a form more palatable to the CASSCF code.  This should */
/*                 be eliminated one day, and the CASSCF code should read    */
/*                 the one-pdm directly.                                     */
/*                                                                           */
/*                 The one-pdm from DETCI is in CI order                     */
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

void onepdm_to_cas(double **onepdm, int *corr_to_pitz, int nso, int npop,
                   int print_lvl, int itap)
{
  PSI_FPTR fptr = 0;            /* index for CAS output file            */
  int cas_opdm_dim = 0;         /* dimension of the onepdm output buff  */
  double *cas_opdm;             /* pointer to in-core onepdm (CI order) */
  int p,q;
  int so_p, so_q, so_pq;

  cas_opdm_dim = (nso * (nso+1))/2;
  cas_opdm = init_array(cas_opdm_dim);

  /* loop through all of the twoel ints and put in Pitzer order in buffer */
  for (p=0; p<npop; p++) {
    so_p = corr_to_pitz[p];

    for (q=0; q<=p; q++) {
      so_q = corr_to_pitz[q];
      so_pq = INDEX(so_p, so_q);
      cas_opdm[so_pq] = onepdm[p][q];
    }
  }


  /*
  ** write information to diskfile and close
  */
  rfile(itap);
  wwritw(itap, (char *) cas_opdm, sizeof(double)*cas_opdm_dim, 
         fptr, &fptr);  
  rclose(itap, 3);

  free(cas_opdm);

}

