/*****************************************************************************/
/* lag_to_cas -    This routine writes out the lagrangian to a               */
/*                 form more palatable to the CASSCF code.                   */
/*                                                                           */
/*                 The lagrangian calculated by this program is in CI order  */
/*                 but CASSCF wants Pitzer order.                            */
/*                 Use pq with p>=q                                          */
/*                                                                           */
/* David Sherrill, April 1998                                                */
/*****************************************************************************/

#include <stdio.h>            /* the standard io library       */
#include <math.h>             /* the C math routines library   */
#include <libciomr/libciomr.h>         /* standard psi routines library */
#include <libqt/qt.h>               /* the quantum trio library      */
#define INDEX(x,y) ((x>y) ? ioff[x] + y : ioff[y] + x)

/*
** declare needed external variables
*/
extern FILE *infile;
extern FILE *outfile;
extern int *ioff;

void lag_to_cas(double **lag, int *corr_to_pitz, int nso,
                int print_lvl, int itap)
{
  PSI_FPTR fptr = 0;            /* index for lagrangian for CAS output file  */
  int cas_lag_dim = 0;          /* dimension of the lagrangian output buff   */
  double **cas_lag;             /* pointer to in-core lagrangian (CI order)  */
  int p,q;
  int so_p, so_q;

  cas_lag_dim = nso * nso;
  cas_lag = block_matrix(nso,nso);

  /* loop through all of the twoel ints and put in Pitzer order in buffer */
  for (p=0; p<nso; p++) {
    so_p = corr_to_pitz[p];

    for (q=0; q<nso; q++) {
      so_q = corr_to_pitz[q];

      cas_lag[so_p][so_q] = lag[p][q];
    }
  }


  /*
  ** write information to diskfile and close
  */
  rfile(itap);
  wwritw(itap, (char *) (cas_lag[0]), sizeof(double)*cas_lag_dim, 
         fptr, &fptr);  
  rclose(itap, 3);

  free_block(cas_lag);

}

