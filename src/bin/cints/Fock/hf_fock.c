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

#include"quartet_data.h"      /* From Default_Ints */
#include"norm_quartet.h"
#include"hash.h"
#include"transmat.h"
#include"int_fjt.h"
#include"schwartz.h"
#include"shell_block_matrix.h"

extern void *hf_fock_thread(void *);
extern pthread_mutex_t fock_mutex;

/*--- To be accessed by all HF Fock threads ---*/
double ****Gskel, ****Gskel_o;       /* Shell-blocked skeleton G matrices */

void hf_fock()
{
  pthread_attr_t thread_attr;
  pthread_t *thread_id;
  
  int count ;
  int dum;
  int g, i, j, k, l, m, ii, jj, kk, ll;
  int si, sj, ni, nj, li, lj, si_g, sj_g;

  double hf_exch = 1.0;                               /* Amount of HF-exchange to be included in G matrix */
  double temp;
  double ****Gfull, ****Gfull_o;       /* Shell-blocked G matrices in AO basis*/
  double ****Gsym, ****Gsym_o;         /* Shell-blocked symmetrized (Gskel + Gskel(transp.)) G matrices */
  double ***ao_type_transmat;

  /*---------------
    Initialization
   ---------------*/
  init_fjt(BasisSet.max_am*4);
  init_libint_base();

  /*------------------------------------------
    Compute integrals for Schwartz inequality
   ------------------------------------------*/
  schwartz_eri();

  /*------------------------------------
    Allocate shell-blocked skeleton G's
   ------------------------------------*/
  Gskel = init_shell_block_matrix();
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    Gskel_o = init_shell_block_matrix();

  thread_id = (pthread_t *) malloc(UserOptions.num_threads*sizeof(pthread_t));
  pthread_attr_init(&thread_attr);
  pthread_attr_setscope(&thread_attr,
			PTHREAD_SCOPE_SYSTEM);
  pthread_mutex_init(&fock_mutex,NULL);
  for(i=0;i<UserOptions.num_threads-1;i++)
    pthread_create(&(thread_id[i]),&thread_attr,
		   hf_fock_thread,(void *)i);
  hf_fock_thread( (void *) (UserOptions.num_threads - 1) );
  for(i=0;i<UserOptions.num_threads-1;i++)
    pthread_join(thread_id[i], NULL);
  free(thread_id);
  
  /*-------------------------------
    Gskel = Gskel + Gskel(transp.)
   -------------------------------*/
  Gsym = init_shell_block_matrix();
  GplusGt(Gskel,Gsym);
  free_shell_block_matrix(Gskel);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    Gsym_o = init_shell_block_matrix();
    GplusGt(Gskel_o,Gsym_o);
    free_shell_block_matrix(Gskel_o);
  }
  
  /*-----------------
    Symmetrize Gskel
   -----------------*/
  if (Symmetry.nirreps > 1) {
    ao_type_transmat = build_transmat(Symmetry.sym_oper, Symmetry.nirreps, BasisSet.max_am);
    Gfull = init_shell_block_matrix();
    for(g=0;g<Symmetry.nirreps;g++) {
      for(si=0;si<BasisSet.num_shells;si++) {
	ni = ioff[BasisSet.shells[si].am];
	li = BasisSet.shells[si].am-1;
	si_g = BasisSet.shells[si].trans_vec[g] - 1;
	for(sj=0;sj<BasisSet.num_shells;sj++) {
	  sj_g = BasisSet.shells[sj].trans_vec[g] - 1;
	  nj = ioff[BasisSet.shells[sj].am];
	  lj = BasisSet.shells[sj].am-1;

	  for(i=0;i<ni;i++)
	    for(j=0;j<nj;j++)
	      Gfull[si_g][sj_g][i][j] += ao_type_transmat[li][g][i]*
					   ao_type_transmat[lj][g][j]*
					   Gsym[si][sj][i][j];
	}
      }
    }
  }
  else
    Gfull = Gsym;

  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    if (Symmetry.nirreps > 1) {
    Gfull_o = init_shell_block_matrix();
    for(g=0;g<Symmetry.nirreps;g++) {
      for(si=0;si<BasisSet.num_shells;si++) {
	ni = ioff[BasisSet.shells[si].am];
	li = BasisSet.shells[si].am-1;
	si_g = BasisSet.shells[si].trans_vec[g] - 1;
	for(sj=0;sj<BasisSet.num_shells;sj++) {
	  sj_g = BasisSet.shells[sj].trans_vec[g] - 1;
	  nj = ioff[BasisSet.shells[sj].am];
	  lj = BasisSet.shells[sj].am-1;

	  for(i=0;i<ni;i++)
	    for(j=0;j<nj;j++)
	      Gfull_o[si_g][sj_g][i][j] += ao_type_transmat[li][g][i]*
					   ao_type_transmat[lj][g][j]*
					   Gsym_o[si][sj][i][j];
	}
      }
    }
    }
    else
      Gfull_o = Gsym_o;
  }
  /*--------------------
    Print out G for now
   --------------------*/
  G = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  shell_block_to_block(Gfull,G);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    Go =  block_matrix(BasisSet.num_ao,BasisSet.num_ao);
    shell_block_to_block(Gfull_o,Go);
  }

  /*---------
    Clean-up
   ---------*/
  free_shell_block_matrix(Gsym);
  if (Symmetry.nirreps > 1) {
    free_shell_block_matrix(Gfull);
  }
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    free_shell_block_matrix(Gsym_o);
    if (Symmetry.nirreps > 1)
      free_shell_block_matrix(Gfull_o);
  }
  free_fjt();

  return;
}

