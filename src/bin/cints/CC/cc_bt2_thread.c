#include<math.h>
#include<stdio.h>
#include<string.h>
#include<memory.h>
#include<stdlib.h>
#include<libqt/qt.h>
#include<libint/libint.h>
#include<pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"quartet_data.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"fjt.h"
#endif


void *cc_bt2_thread(void *tnum_ptr)
{
  int thread_num = (int) tnum_ptr;
  const double toler = UserOptions.cutoff;
  const double m_sqrt1_2 = 1/sqrt(2.0);

  /*--- Various data structures ---*/
  struct shell_pair *sp_ij, *sp_kl;
  struct unique_shell_pair *usp_ij,*usp_kl;
  Libint_t Libint;
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;
#endif

  int ij, kl, ik, jl, ijkl;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int i, m, p, q, r, s;
  int p_abs, q_abs, r_abs, s_abs;
  int si, sj, sk, sl ;
  int sii, sjj, skk, sll , slll;
  int num_ij, swap_ij_kl;
  int pi, pj, pk, pl;
  int max_pj, max_pl;
  int *sj_arr, *sk_arr, *sl_arr;
  int sr_fao, ss_fao, sp_fao, sq_fao;
  int usii,usjj,uskk,usll,usi,usj,usk,usl,usij;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int lambda_T = 1;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

  int size;
  int max_class_size;
  int max_cart_class_size;

  int np_i, np_j, np_k, np_l;
  int nr, ns, np, nq;
  int num_prim_comb;

  int num_ibatch, num_i_per_ibatch, ibatch, ibatch_length;
  int imin, imax, jmin;
  int max_bf_per_shell;
  int mo_i, mo_j, mo_a, mo_b, mo_ij;
  int ia;
  int rs_offset, rsi_offset, rsp_offset;

  double AB2, CD2;

  double *raw_data;             /* pointer to the unnormalized taregt quartet of integrals */
  double *data;                 /* pointer to the transformed normalized target quartet of integrals */
#ifdef NONDOUBLE_INTS
  REALTYPE *target_ints;            /* Pointer to the location of the target quartet on the stack of
			 	    integrals quartets if libint.a is using other than regular doubles */
#endif

  double *rspq_ptr;
  double temp;
  double *mo_vec;
  double *rsiq_buf;             /* buffer for (rs|iq) integrals, where r,s run over shell sets,
				   i runs over I-batch, q runs over all AOs */
  double *rsi_row, *i_row;
  double *ia_block_ptr;
  double *rsia_buf;             /* buffer for (rs|ia) integrals, where r,s run over shell sets,
				   i runs over I-batch, q runs over all AOs */
  double *jsi_row;
  double *jbi_row;
  double iajb, ibja, pfac, k0ijab, k1ijab, eijab, e0, e1;

  double temp1,temp2,*iq_row,*ip_row;
  int rs,qrs;

  double *scratch_buf;          /* scratch used in permuting bra and ket */

  /*---------------
    Initialization
   ---------------*/
#ifndef USE_TAYLOR_FM
  init_fjt_table(&fjt_table);
#endif
  
  /*-------------------------
    Allocate data structures
   -------------------------*/
  max_bf_per_shell = ioff[BasisSet.max_am];
  max_cart_class_size = (max_bf_per_shell)*
                        (max_bf_per_shell)*
                        (max_bf_per_shell)*
                        (max_bf_per_shell);
  max_num_unique_quartets = Symmetry.max_stab_index*
                            Symmetry.max_stab_index*
                            Symmetry.max_stab_index;
  sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  if (Symmetry.nirreps > 1)
    max_class_size = max_cart_class_size;
  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  init_libint(&Libint, BasisSet.max_am-1, max_num_prim_comb);

#ifdef NONDOUBLE_INTS
  raw_data = init_array(max_cart_class_size);
#endif

  
/*-----------------------------------
  generate all unique shell quartets
 -----------------------------------*/
  /*--- "unique" R,S loop ---*/
  usij = 0;
  for (usii=0; usii<Symmetry.num_unique_shells; usii++)
    for (usjj=0; usjj<=usii; usjj++, usij++) {
      /*--- Decide if this thread will do this ---*/
      if ( usij%UserOptions.num_threads != thread_num )
	continue;
      usi = usii; usj = usjj;
      /*--- As usual, swap order usi and usj according to their angular momenta ---*/
      if(BasisSet.shells[Symmetry.us2s[usi]].am < BasisSet.shells[Symmetry.us2s[usj]].am){
	dum = usi;
	usi = usj;
	usj = dum;
      }
	
      sii = Symmetry.us2s[usi];
      sjj = Symmetry.us2s[usj];
      if (Symmetry.nirreps > 1) {
	stab_i = Symmetry.atom_positions[BasisSet.shells[sii].center-1];
	stab_j = Symmetry.atom_positions[BasisSet.shells[sjj].center-1];
	stab_ij = Symmetry.GnG[stab_i][stab_j];
	R_list = Symmetry.dcr[stab_i][stab_j];
	num_ij = Symmetry.dcr_dim[stab_i][stab_j];
      }
      else
	num_ij = 1;
      
      /*--- R,S loop ---*/
      for(dcr_ij=0;dcr_ij<num_ij;dcr_ij++) {
	if (Symmetry.nirreps > 1)
	  R = R_list[dcr_ij];
	else
	  R = 0;
	si = sii;
	sj = BasisSet.shells[sjj].trans_vec[R]-1;
	
	/*--- "Unique" R,S loop ---*/
	for (uskk=0; uskk<Symmetry.num_unique_shells; uskk++)
	  for (usll=0; usll<=uskk; usll++){
	    
	    /*--- For each combination of unique shells generate "petit list" of shells ---*/
	    usk = uskk; usl = usll;
	    /*--- As usual, swap order usk and usl according to their angular momenta ---*/
	    if(BasisSet.shells[Symmetry.us2s[usk]].am < BasisSet.shells[Symmetry.us2s[usl]].am){
	      dum = usk;
	      usk = usl;
	      usl = dum;
	    }
	    /*--- DO NOT SWAP bra and ket at this time. Do it later, in the main loop ---*/
	    if(BasisSet.shells[Symmetry.us2s[usi]].am + BasisSet.shells[Symmetry.us2s[usj]].am >
	       BasisSet.shells[Symmetry.us2s[usk]].am + BasisSet.shells[Symmetry.us2s[usl]].am)
	      swap_ij_kl = 1;
	    else
	      swap_ij_kl = 0;

	    skk = Symmetry.us2s[usk];
	    sll = Symmetry.us2s[usl];
	    if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
	      /*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
	      stab_k = Symmetry.atom_positions[BasisSet.shells[skk].center-1];
	      stab_l = Symmetry.atom_positions[BasisSet.shells[sll].center-1];
	      stab_kl = Symmetry.GnG[stab_k][stab_l];
	      S_list = Symmetry.dcr[stab_k][stab_l];
	      T_list = Symmetry.dcr[stab_ij][stab_kl];
	      lambda_T = Symmetry.nirreps/Symmetry.dcr_deg[stab_ij][stab_kl];

	      memset(sj_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sk_arr,0,sizeof(int)*max_num_unique_quartets);
	      memset(sl_arr,0,sizeof(int)*max_num_unique_quartets);
	      count = 0;
		
	      for(dcr_ijkl=0;dcr_ijkl<Symmetry.dcr_dim[stab_ij][stab_kl];dcr_ijkl++){
		T = T_list[dcr_ijkl];
		sk = BasisSet.shells[skk].trans_vec[T]-1;
		slll = BasisSet.shells[sll].trans_vec[T]-1;
		for(dcr_kl=0;dcr_kl<Symmetry.dcr_dim[stab_k][stab_l];dcr_kl++) {
		  S = S_list[dcr_kl];
		  sl = BasisSet.shells[slll].trans_vec[S]-1;
		  
		  total_am = BasisSet.shells[si].am +
		    BasisSet.shells[sj].am +
		    BasisSet.shells[sk].am +
		    BasisSet.shells[sl].am;
		  /*-------------------------------------------------------------
		    Obviously redundant or zero cases should be eliminated here!
		    Right now only zero case is eliminated. Redundancies arising
		    in DCD approach when usi == usj etc. may be eliminated too
		    but lambda_T will have to be replaced by an array (it won't
		    the same for every shell quartet in petite list anymore).
		   -------------------------------------------------------------*/
		  if(!(total_am%2)||
		     (BasisSet.shells[si].center!=BasisSet.shells[sj].center)||
		     (BasisSet.shells[sj].center!=BasisSet.shells[sk].center)||
		     (BasisSet.shells[sk].center!=BasisSet.shells[sl].center)) {
		    sj_arr[count] = sj;
		    sk_arr[count] = sk;
		    sl_arr[count] = sl;
		    count++;
		  }
		}
	      } /* petite list is ready to be used */
	      num_unique_quartets = count;
	    }
	    else { /*--- C1 symmetry case ---*/
	      total_am = BasisSet.shells[si].am +
		BasisSet.shells[usj].am +
		BasisSet.shells[usk].am +
		BasisSet.shells[usl].am;
	      if(!(total_am%2)||
		 (BasisSet.shells[si].center!=BasisSet.shells[usj].center)||
		 (BasisSet.shells[usj].center!=BasisSet.shells[usk].center)||
		 (BasisSet.shells[usk].center!=BasisSet.shells[usl].center)) {
		num_unique_quartets = 1;
		sj_arr[0] = usj;
		sk_arr[0] = usk;
		sl_arr[0] = usl;
	      }
	      else
		num_unique_quartets = 0;
	    }
	    
	    /*----------------------------------
	      Compute the nonredundant quartets
	     ----------------------------------*/
	    for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
	      si = sii;
	      sj = sj_arr[plquartet];
	      sk = sk_arr[plquartet];
	      sl = sl_arr[plquartet];
	      /*--- As usual, we have to order bra-ket so that ket has the largest angular momentum */
	      if (swap_ij_kl) {
		dum = si;
		si = sk;
		sk = dum;
		dum = sj;
		sj = sl;
		sl = dum;
	      }
	      np_i = BasisSet.shells[si].n_prims;
	      np_j = BasisSet.shells[sj].n_prims;
	      np_k = BasisSet.shells[sk].n_prims;
	      np_l = BasisSet.shells[sl].n_prims;
	      orig_am[0] = BasisSet.shells[si].am-1;
	      orig_am[1] = BasisSet.shells[sj].am-1;
	      orig_am[2] = BasisSet.shells[sk].am-1;
	      orig_am[3] = BasisSet.shells[sl].am-1;
	      am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];
	      
	      sp_ij = &(BasisSet.shell_pairs[si][sj]);
	      sp_kl = &(BasisSet.shell_pairs[sk][sl]);
	      
	      Libint.AB[0] = sp_ij->AB[0];
	      Libint.AB[1] = sp_ij->AB[1];
	      Libint.AB[2] = sp_ij->AB[2];
	      Libint.CD[0] = sp_kl->AB[0];
	      Libint.CD[1] = sp_kl->AB[1];
	      Libint.CD[2] = sp_kl->AB[2];
	      
	      AB2 = Libint.AB[0]*Libint.AB[0]+Libint.AB[1]*Libint.AB[1]+Libint.AB[2]*Libint.AB[2];
	      CD2 = Libint.CD[0]*Libint.CD[0]+Libint.CD[1]*Libint.CD[1]+Libint.CD[2]*Libint.CD[2];
	      
	      /*--- Compute data for primitive quartets here ---*/
	      num_prim_comb = 0;
	      for (pi = 0; pi < np_i; pi++) {
		max_pj = (si == sj) ? pi+1 : np_j;
		for (pj = 0; pj < max_pj; pj++) {
		  m = (1 + (si == sj && pi != pj));
		  for (pk = 0; pk < np_k; pk++) {
		    max_pl = (sk == sl) ? pk+1 : np_l;
		    for (pl = 0; pl < max_pl; pl++){
		      n = m * (1 + (sk == sl && pk != pl));
#ifdef USE_TAYLOR_FM
		      quartet_data(&(Libint.PrimQuartet[num_prim_comb++]), NULL, AB2, CD2,
				   sp_ij, sp_kl, am, pi, pj, pk, pl, n*lambda_T);
#else
		      quartet_data(&(Libint.PrimQuartet[num_prim_comb++]), &fjt_table, AB2, CD2,
				   sp_ij, sp_kl, am, pi, pj, pk, pl, n*lambda_T);
#endif
		    }
		  }
		}
	      }

	      /*--- Compute the integrals ---*/
	      if (am) {
#ifdef NONDOUBLE_INTS
		size = ioff[BasisSet.shells[si].am]*ioff[BasisSet.shells[sj].am]*
		  ioff[BasisSet.shells[sk].am]*ioff[BasisSet.shells[sl].am];
		target_ints = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint, num_prim_comb);
		for(i=0;i<size;i++)
		  raw_data[i] = (double) target_ints[i];
#else
		raw_data = build_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libint, num_prim_comb);
#endif
		/* No need to transforms integrals to sph. harm. basis */
		data = norm_quartet(raw_data, NULL, orig_am, 0);
	      }
	      else {
		temp = 0.0;
		for(p=0;p<num_prim_comb;p++)
		  temp += (double) Libint.PrimQuartet[p].F[0];
#ifdef NONDOUBLE_INTS
		raw_data[0] = temp;
		data = raw_data;
#else
		Libint.int_stack[0] = temp;
		data = Libint.int_stack;
#endif
	      }

	      /*--- swap bra and ket back to the original order if needed ---*/
	      if (swap_ij_kl) {
		dum = si;
		si = sk;
		sk = dum;
		dum = sj;
		sj = sl;
		sl = dum;
		sr_fao = BasisSet.shells[si].fao - 1;
		ss_fao = BasisSet.shells[sj].fao - 1;
		sp_fao = BasisSet.shells[sk].fao - 1;
		sq_fao = BasisSet.shells[sl].fao - 1;
		nr = ioff[BasisSet.shells[si].am];
		ns = ioff[BasisSet.shells[sj].am];
		np = ioff[BasisSet.shells[sk].am];
		nq = ioff[BasisSet.shells[sl].am];
		if (am /*&& (orig_am[0] + orig_am[1] != 0)*/) {
		  /*		    timer_on("Pre1Swap");*/
		  /*--- (pq|rs) -> (rs|pq) ---*/
		  ijkl_to_klij(data,scratch_buf,np*nq,nr*ns);
		  data = scratch_buf;
		  /*		    timer_off("Pre1Swap");*/
		}
	      }
	      else {
		sr_fao = BasisSet.shells[si].fao - 1;
		ss_fao = BasisSet.shells[sj].fao - 1;
		sp_fao = BasisSet.shells[sk].fao - 1;
		sq_fao = BasisSet.shells[sl].fao - 1;
		nr = ioff[BasisSet.shells[si].am];
		ns = ioff[BasisSet.shells[sj].am];
		np = ioff[BasisSet.shells[sk].am];
		nq = ioff[BasisSet.shells[sl].am];
	      }

	    } /* end of RSPQ loop */
	  } /* end of "unique" P,S loop */
      } /* end of R,S loop */
    } /* end of "unique" R,S loop */

  /*---------
    Clean-up
   ---------*/
#ifdef NONDOUBLE_INTS
  free(raw_data);
#endif
  free_libint(&Libint);
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
#ifndef USE_TAYLOR_FM
  free_fjt_table(&fjt_table);
#endif


  return;
}


