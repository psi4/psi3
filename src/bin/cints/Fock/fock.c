#include<math.h>
#include<stdio.h>
#include<string.h>
#include<memory.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<psio.h>
#include<libint.h>
#include<pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"

#include"read_scf_opdm.h"
#include"read_scf_occ_evec.h"
#include"shell_block_matrix.h"
#include"hf_fock.h"
#include"xc_fock.h"

pthread_mutex_t fock_mutex;            /* Lock on the global AO matrix */

void fock()
{
  int dum;
  int n, num;
  int i, j, k, l;

  int nstri;
  double temp;
  double **tmpmat1;
  double *Gtri, *Gtri_o;  /* Total and open-shell G matrices 
			     and lower triagonal form in SO basis */

  /*----------------------------------------
    Read in the difference HF/DFT densities
   ----------------------------------------*/
  read_scf_opdm();

  /*-------------------------------------------
    Compute HF contribution to the Fock matrix
   -------------------------------------------*/
  hf_fock();
  
  /*-----------------------------------
    Do numerical interation for KS DFT
   -----------------------------------*/
  if(UserOptions.make_dft){
    /*--- Read in the SCF eigenvector density ---*/
      read_scf_occ_evec();
    /*-- Compute exch+corr contribution to the Fock matrix ---*/
    xc_fock();
  }

  return;
}

