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
#include"read_scf_opdm.h"
#include"int_fjt.h"
#include"schwartz.h"
#include"shell_block_matrix.h"

extern void *fock_thread(void *);
pthread_mutex_t fock_mutex;

/*--- To be accessed by all threads ---*/
double ****Gskel, ****Gskel_o;       /* Shell-blocked skeleton G matrices */

void fock()
{
  pthread_attr_t thread_attr;
  pthread_t *thread_id;
  
  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int g, i, j, k, l, m, ii, jj, kk, ll;
  int a, b, c, d;
  int si;                                  /* GCC compiler screwes up if static is left out */
  int sj, sk, sl, si_g, sj_g;
  int sii, sjj, skk, sll , slll;
  int sij, skl, sijkl;
  int pi, pj, pk, pl ;
  int max_pj, max_pl;
  int pii, pjj, pkk, pll;
  int upk, num_unique_pk;
  int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
  int *si_arr, *sj_arr, *sk_arr, *sl_arr, *key_arr;
  int usii,usjj,uskk,usll,usi,usj,usk,usl;
  int usi_eq_usj, usi_eq_usk, usi_eq_usl, usj_eq_usl, usk_eq_usj, usk_eq_usl;
  int usij_eq_uskl, usik_eq_usjl, usil_eq_uskj;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

  int class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, li, lj;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;
  int key, key1, key2, key3, new_quartet, htable_ptr;
  int nstri;

  double hf_exch = 1.0;                               /* Amount of HF-exchange to be included in G matrix */
  double AB2, CD2;
  double *data;
  double temp;
  double **tmpmat1;
  double *Gtri, *Gtri_o;               /* Total and open-shell G matrices and lower triagonal form
					  in SO basis */
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

  /*-----------------------------
    Read in the HF/DFT densities
   -----------------------------*/
  read_scf_opdm();
  
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
		   fock_thread,(void *)i);
  fock_thread( (void *) (UserOptions.num_threads - 1) );
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

  /*------------------
     DFT Procedure
    ------------------*/
  
  
  /*--------------------
    Print out G for now
   --------------------*/
  G = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  shell_block_to_block(Gfull,G);
/*  fprintf(outfile,"  Closed-shell Fock matrix in AO basis:\n");
  print_mat(G,BasisSet.num_ao,BasisSet.num_ao,outfile);*/
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    Go =  block_matrix(BasisSet.num_ao,BasisSet.num_ao);
    shell_block_to_block(Gfull_o,Go);
/*    fprintf(outfile,"  Open-shell Fock matrix in AO basis:\n");
    print_mat(Go,BasisSet.num_ao,BasisSet.num_ao,outfile);*/
  }
  if (Symmetry.nirreps > 1 || BasisSet.puream) {
    tmpmat1 = block_matrix(Symmetry.num_so,BasisSet.num_ao);
    mmult(Symmetry.usotao,0,G,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
    mmult(tmpmat1,0,Symmetry.usotao,1,G,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
    if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      mmult(Symmetry.usotao,0,Go,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,Go,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
    }
    free_block(tmpmat1);
  }
/*  fprintf(outfile,"  Closed-shell Fock matrix in SO basis:\n");
  print_mat(G,Symmetry.num_so,Symmetry.num_so,outfile);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    fprintf(outfile,"  Open-shell Fock matrix in SO basis:\n");
    print_mat(Go,Symmetry.num_so,Symmetry.num_so,outfile);
  }*/

  /*-------------------------
    Write G-matrices to disk
   -------------------------*/
  nstri = ioff[Symmetry.num_so];
  Gtri = init_array(nstri);
  sq_to_tri(G,Gtri,Symmetry.num_so);
  free_block(G);
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_NEW);
  switch (UserOptions.reftype) {
  case rohf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      psio_write_entry(IOUnits.itapDSCF, "Open-shell G-matrix", (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri_o);
  case rhf:
      psio_write_entry(IOUnits.itapDSCF, "Total G-matrix", (char *) Gtri, sizeof(double)*nstri);
      free(Gtri);
      break;

  case uhf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      /*--- Form alpha and beta Fock matrices first and then write them out ---*/
      for(i=0;i<nstri;i++) {
	temp = Gtri[i] + Gtri_o[i];
	Gtri[i] = Gtri[i] - Gtri_o[i];
	Gtri_o[i] = temp;
      }
      psio_write_entry(IOUnits.itapDSCF, "Alpha G-matrix", (char *) Gtri, sizeof(double)*nstri);
      psio_write_entry(IOUnits.itapDSCF, "Beta G-matrix", (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri);
      free(Gtri_o);
      break;
  }
  psio_close(IOUnits.itapDSCF, 1);

  
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

