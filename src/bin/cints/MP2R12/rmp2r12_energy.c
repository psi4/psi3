#include<math.h>
#include<stdio.h>
#include<string.h>
#include<memory.h>
#include<stdlib.h>
#include<libciomr.h>
#include<qt.h>
#include<iwl.h>
#include<libint.h>
#include<libr12.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"r12_quartet_data.h"
#include"norm_quartet.h"
#include"int_fjt.h"
#include"quartet_permutations.h"


/*-------------------------------
  Explicit function declarations
 -------------------------------*/
void make_transqt_arrays(int **first, int **last, int **fstocc, int **lstocc, int **occ, int **act2fullQTS, int **ioff3);

/*-------------------------------------------------------
  Algorithm

  Loop over I batches (batch size num_i_per_batch) of active DOCC

    Loop over all symmetry-unique shells UR, US<=UR
      Find all symmetry-distinct shell doublets resulting from (UR US|
      Loop over all resulting shell doublets R, S

        Loop over all symmetry-unique shells UP, UQ<=UP
	  Find all symmetry-distinct shell quartets resulting from (R S|UP UQ)
	  Loop over the resulting set of P, Q doublets
            Evaluate (RS|PQ), (RS|r12|PQ), and (RS|[r12,T1]|PQ)
	    Loop over p in P, q in Q, r in R, s in S, i in I
	      (rs|iq) += Cpi * (rs|pq)
	      (rs|ip) += Cqi * (rs|pq)
	      same for (rs|r12|pq) and (rs|[r12,T1]|pq)
	    End p, q, r, s, i loop
          End P,Q loop
        End UP, UQ loop

        Loop over r in R, s in S
          Loop over q < nao, x < num_mo, i in I
	    (rs|ix) += Cxs * (rs|is)
	    same for (rs|r12|is) and (rs|[r12,T1]|is)
	  End q, x, i
	  Loop over i in I, x < num_mo, j <= i
	    (js|ix) += Cjr * (rs|ix)
	    (jr|ix) += Cjs * (rs|ix)
	    same for (rs|r12|ix), but
	    (js|[r12,T1]|ix) += Cjr * (rs|[r12,T1]|ix)
	    (jr|[r12,T1]|ix) -= Cjs * (rs|[r12,T1]|ix)                   <---- Note the minus sign here!!!
	  End i, x, j loop
        End r, s loop

      End R, S loop
    End UR, US loop

    Loop over i in I, j <= i
      Loop over r < nao, x < num_mo, y < num_mo
        (jy|ix) += Cys * (js|ix)
	same for (js|r12|ix) and (js|[r12,T1]|ix)
      End r, x, y loop
    End i, j loop

  End I loop
  
 -------------------------------------------------------*/


void rmp2r12_energy()
{
  const double toler = UserOptions.cutoff;

  struct iwlbuf MOBuf[3];
  struct shell_pair *sp_ij, *sp_kl;
  Libr12_t Libr12;
  double_array_t fjt_table;

  static char *te_operator[] = { "1/r12", "r12", "[r12,T1]" };
  int total_te_count[NUM_TE_TYPES] = {0, 0, 0, 0};
  int ij, kl, ik, jl, ijkl;
  int count ;
  int dum;
  int te_type;
  int n, num[NUM_TE_TYPES];
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  int i, j, m, p, q, r, s, x, y;
  int isym, jsym, xsym, ysym;
  int p_abs, q_abs, r_abs, s_abs;
  int si, sj, sk, sl;
  int sii, sjj, skk, sll , slll;
  int num_ij, swap_ij_kl;
  int pi, pj, pk, pl ;
  int *sj_arr, *sk_arr, *sl_arr;
  int sr_fao, ss_fao, sp_fao, sq_fao;
  int usii,usjj,uskk,usll,usi,usj,usk,usl;
  int stab_i,stab_j,stab_k,stab_l,stab_ij,stab_kl;
  int *R_list, *S_list, *T_list;
  int R,S,T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int lambda_T = 1;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

  int size, class_size;
  int max_class_size;
  int max_cart_class_size;

  int np_i, np_j, np_k, np_l;
  int nr, ns, np, nq;

  int num_ibatch, num_i_per_ibatch, ibatch, ibatch_length;
  int imin, imax, jmin;
  int max_bf_per_shell;
  int mo_i, mo_j, mo_x, mo_y;
  int rs, qrs;
  int rs_offset, rsi_offset, rsp_offset;
  int num_prim_comb;

  double Emp2 = 0.0;
  double Emp2r12 = 0.0;
  double AB2, CD2;
  double *raw_data[NUM_TE_TYPES];             /* pointers to the unnormalized taregt quartet of integrals */
  double *data[NUM_TE_TYPES];                 /* pointers to the transformed normalized target quartet of integrals */

  double temp;
  double ssss, ss_r12_ss;
  double *rspq_ptr;
  double *mo_vec;
  double *rsiq_buf[NUM_TE_TYPES];             /* buffer for (rs|iq) integrals, where r,s run over shell sets,
						 i runs over I-batch, q runs over all AOs */
  double *rsi_row;
  double **ix_buf;                            /* buffer for one |ix) ket */
  double *i_row;
  double **jy_buf;                            /* buffer for one |jy) ket */
  double *jsix_buf[NUM_TE_TYPES];             /* buffer for (js|ia) integrals, where j runs over all d.-o. MOs,
						 s runs over all AOs, i - over I-batch, x - over all MOs */
  double *jsi_row;
  double *jyix_buf[NUM_TE_TYPES];             /* buffer contains all MP2-R12/A-type integrals */
  double *jyi_row;
  double *ixjy_buf[NUM_TE_TYPES];                 /* just pointers to (ix|jy) integrals */
  double temp1,temp2,*iq_row,*ip_row;

  int *first, *last;                          /* first and last absolute (Pitzer) orbital indices in symblk */
  int *fstocc, *lstocc;                       /* first and last occupied indices in Pitzer ordering for each symblk */
  int *occ;                                   /* Pitzer to "full"(no frozen core) QTS index mapping */
  int *act2fullQTS;                           /* Maps "active"(taking frozen core into account) QTS into "frozen" QTS index */
  int *ioff3;                                 /* returns pointers to the beginning of rows in a rectangular matrix */

  /*---------------
    Initialization
   ---------------*/
  init_fjt(BasisSet.max_am*4+1);
  init_libr12_base();
  iwl_buf_init(&MOBuf[0], IOUnits.itapERI_MO, toler, 0, 0);
  iwl_buf_init(&MOBuf[1], IOUnits.itapR12_MO, toler, 0, 0);
  iwl_buf_init(&MOBuf[2], IOUnits.itapR12T2_MO, toler, 0, 0);
  make_transqt_arrays(&first, &last, &fstocc, &lstocc, &occ, &act2fullQTS, &ioff3);
  timer_init();
  fprintf(outfile,"  Performing AO->MO integral tranformation for RHF MP2-R12/A via direct algorithm\n");
  
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
#ifdef NONDOUBLE_INTS
  for(i=0;i<NUM_TE_TYPES;i++)
    raw_data[i] = init_array(max_cart_class_size);
#endif
  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  UserOptions.memory -= init_libr12(&Libr12,max_num_prim_comb);
  init_fjt_table(&fjt_table);


  /*---
    Minimum number of I-batches - 
    take sizes of rsiq_buf, jsix_buf, and
    jyix_buf into account
   ---*/
  num_i_per_ibatch = UserOptions.memory / (NUM_TE_TYPES*
					   (BasisSet.num_ao*max_bf_per_shell*max_bf_per_shell +
					    MOInfo.num_mo*
					    (MOInfo.nactdocc*BasisSet.num_ao +
					     MOInfo.nactdocc*MOInfo.num_mo)));
  if (num_i_per_ibatch > MOInfo.nactdocc)
    num_i_per_ibatch = MOInfo.nactdocc;
  if (num_i_per_ibatch < 1)
    punt("Not enough memory for direct MP2-R12/A transformation");
  num_ibatch = (MOInfo.nactdocc + num_i_per_ibatch - 1) / num_i_per_ibatch;
  /*--- Recompute number of MOs per I-batch ---*/
  num_i_per_ibatch = (MOInfo.nactdocc + num_ibatch - 1) / num_ibatch;
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    rsiq_buf[te_type] = init_array(num_i_per_ibatch*BasisSet.num_ao*
				   max_bf_per_shell*max_bf_per_shell);
    jsix_buf[te_type] = init_array(MOInfo.nactdocc*BasisSet.num_ao*
				   num_i_per_ibatch*MOInfo.num_mo);
    jyix_buf[te_type] = init_array(MOInfo.nactdocc*MOInfo.num_mo*
				   num_i_per_ibatch*MOInfo.num_mo);
  }
  ix_buf = block_matrix(num_i_per_ibatch,MOInfo.num_mo);
  jy_buf = block_matrix(MOInfo.ndocc,MOInfo.num_mo);
  fprintf(outfile,"  Using %d %s\n\n",num_ibatch, (num_ibatch == 1) ? "pass" : "passes");


/*-----------------------------------
  generate all unique shell quartets
 -----------------------------------*/
  /*--- I-batch loop ---*/
  for (ibatch=0;ibatch<num_ibatch;ibatch++) {
    imin = ibatch * num_i_per_ibatch + MOInfo.nfrdocc;
    imax = MIN( imin+num_i_per_ibatch , MOInfo.ndocc );
    ibatch_length = imax - imin;
    jmin = MOInfo.nfrdocc;
    fprintf(outfile,"  Pass #%d, MO %d through MO %d\n",ibatch,imin+1,imax);
    fflush(outfile);

    /*--- "unique" R,S loop ---*/
    for (usii=0; usii<Symmetry.num_unique_shells; usii++)
      for (usjj=0; usjj<=usii; usjj++) {
	usi = usii; usj = usjj;

	/*----------------------------------------------------------------------
	  NOTE on swapping usi, usj, etc.:

	  For 2-electron integrals of Hermitian operators it does not matter
	  if we swap si and sj, or sk and sl, or even si,sj and sk,sl. It
	  matters here though since [r12,T1] is non-Hermitian! What we want
	  in the end are the integrals or [r12,T1] operator of this type:
	  (si sj|[r12,T1]|sk sl). If we have to swap bra and ket in the process
	  we will end up with having to compute (sk sl|[r12,T2]|si sj) instead.
	  Therefore if we have to swap bra and ket (swap_ij_kl = 1), then we
	  have to take [r12,T2] integral instead and swap it's bra and ket back.
	 ----------------------------------------------------------------------*/
	
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

	  /*--- "Unique" P,Q loop ---*/
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

		Libr12.ShellQuartet.AB[0] = sp_ij->AB[0];
		Libr12.ShellQuartet.AB[1] = sp_ij->AB[1];
		Libr12.ShellQuartet.AB[2] = sp_ij->AB[2];
		Libr12.ShellQuartet.CD[0] = sp_kl->AB[0];
		Libr12.ShellQuartet.CD[1] = sp_kl->AB[1];
		Libr12.ShellQuartet.CD[2] = sp_kl->AB[2];
		Libr12.ShellQuartet.AC[0] = Molecule.centers[BasisSet.shells[si].center-1].x-
			Molecule.centers[BasisSet.shells[sk].center-1].x;
		Libr12.ShellQuartet.AC[1] = Molecule.centers[BasisSet.shells[si].center-1].y-
			Molecule.centers[BasisSet.shells[sk].center-1].y;
		Libr12.ShellQuartet.AC[2] = Molecule.centers[BasisSet.shells[si].center-1].z-
			Molecule.centers[BasisSet.shells[sk].center-1].z;
		Libr12.ShellQuartet.ABdotAC = Libr12.ShellQuartet.AB[0]*Libr12.ShellQuartet.AC[0]+
					      Libr12.ShellQuartet.AB[1]*Libr12.ShellQuartet.AC[1]+
					      Libr12.ShellQuartet.AB[2]*Libr12.ShellQuartet.AC[2];
		Libr12.ShellQuartet.CDdotCA = -1.0*(Libr12.ShellQuartet.CD[0]*Libr12.ShellQuartet.AC[0]+
						    Libr12.ShellQuartet.CD[1]*Libr12.ShellQuartet.AC[1]+
						    Libr12.ShellQuartet.CD[2]*Libr12.ShellQuartet.AC[2]);
		AB2 = Libr12.ShellQuartet.AB[0]*Libr12.ShellQuartet.AB[0]+
		      Libr12.ShellQuartet.AB[1]*Libr12.ShellQuartet.AB[1]+
		      Libr12.ShellQuartet.AB[2]*Libr12.ShellQuartet.AB[2];
		CD2 = Libr12.ShellQuartet.CD[0]*Libr12.ShellQuartet.CD[0]+
		      Libr12.ShellQuartet.CD[1]*Libr12.ShellQuartet.CD[1]+
		      Libr12.ShellQuartet.CD[2]*Libr12.ShellQuartet.CD[2];

		/*--------------------------------
		  contract by primitives out here
		 --------------------------------*/
		num_prim_comb = 0;
		for (pi = 0; pi < np_i; pi++)
		  for (pj = 0; pj < np_j; pj++)
		    for (pk = 0; pk < np_k; pk++)
		      for (pl = 0; pl < np_l; pl++){
			r12_quartet_data(&(Libr12.PrimQuartet[num_prim_comb++]), &fjt_table, AB2, CD2,
					 sp_ij, sp_kl, am, pi, pj, pk, pl, lambda_T);
		      }

		if (am) {
		  build_r12_grt[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libr12,num_prim_comb);
		  if (swap_ij_kl)
		    /*--- (usi usj|[r12,T1]|usk usl) = (usk usl|[r12,T2]|usi usj) ---*/
		    Libr12.te_ptr[2] = Libr12.te_ptr[3];
#ifdef NONDOUBLE_INTS
		  size = ioff[BasisSet.shells[si].am]*ioff[BasisSet.shells[sj].am]*
			 ioff[BasisSet.shells[sk].am]*ioff[BasisSet.shells[sl].am];
		  for(i=0;i<NUM_TE_TYPES-1;i++)
		    for(j=0;j<size;j++)
		      raw_data[i][j] = (double) Libr12.te_ptr[i][j];
#else
		  for(i=0;i<NUM_TE_TYPES-1;i++)
		    raw_data[i] = Libr12.te_ptr[i];
#endif
		  /*--- Just normalize the integrals ---*/
		  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++)
		    data[te_type] = norm_quartet(raw_data[te_type], NULL, orig_am, 0);
		}
		else {
		  ssss = 0.0;
		  ss_r12_ss = 0.0;
		  for(p=0;p<num_prim_comb;p++) {
		    ssss += (double) Libr12.PrimQuartet[p].F[0];
		    ss_r12_ss += (double) Libr12.PrimQuartet[p].ss_r12_ss;
		  }
		  build_r12_grt[0][0][0][0](&Libr12,num_prim_comb);
#ifdef NONDOUBLE_INTS
		  raw_data[0][0] = ssss;
		  raw_data[1][0] = ss_r12_ss;
		  raw_data[2][0] = (double) Libr12.te_ptr[2][0];
		  raw_data[3][0] = (double) Libr12.te_ptr[3][0];
		  data[0] = raw_data[0];
		  data[1] = raw_data[1];
		  data[2] = raw_data[2];
		  data[3] = raw_data[3];
#else
		  Libr12.int_stack[2] = Libr12.te_ptr[2][0];
		  Libr12.int_stack[3] = Libr12.te_ptr[3][0];
		  Libr12.int_stack[0] = ssss;
		  Libr12.int_stack[1] = ss_r12_ss;
		  data[0] = Libr12.int_stack;
		  data[1] = Libr12.int_stack+1;
		  data[2] = Libr12.int_stack+2;
		  data[3] = Libr12.int_stack+3;
#endif
		}

		/*---
		  Swap bra and ket back to the reversed order (PQ|RS) if not done yet
		  Need this to make Step 1 a (fast) vector operation!
		  NOTE: This changes only the way integrals are stored! No need to
		  worry about non-hermiticity of [r12,T1] here!!!
		  ---*/
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
		  /*---
		    No need to swap bra and ket since the quartet
		    was already computed in reverse order
		   ---*/
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
		  /*---
		    Need to swap bra and ket, but do
		    the actual permutation in te_type-loop
		   ---*/
		}

		/*--- step 1 of the transformation ---*/
		timer_on("Step 1");
		for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
		  /*---
		    if bra and ket were not swapped before computing
		    the quartet, swap them now so that Step 1 is
		    a vector operation
		   ---*/
		  if ((!swap_ij_kl) && am) {
		    timer_on("Pre1Swap");
		    /*--- (rs|pq) -> (pq|rs) ---*/
		    data[te_type] = ijkl_to_klij(data[te_type],nr*ns,np*nq);
		    timer_off("Pre1Swap");
		  }
		  if (usk != usl)
		    for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		      mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
		      rspq_ptr = data[te_type];
		      for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			for(q=0,q_abs=sq_fao;
			    q<nq;
			    q++,q_abs++) {
			  ip_row = rsiq_buf[te_type] + ((mo_i*BasisSet.num_ao + p_abs) * nr * ns);
			  temp1 = mo_vec[q_abs];
#if !USE_BLAS
			  for(rs=0;rs<nr*ns;rs++,rspq_ptr++,iq_row++,ip_row++) {
			    (*ip_row) += temp1 * (*rspq_ptr);
			  }
#else
			  C_DAXPY(nr*ns,temp1,rspq_ptr,1,ip_row,1);
			  rspq_ptr += nr*ns;
#endif
			}
		      }
		      rspq_ptr = data[te_type];
		      for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			temp = mo_vec[p_abs];
			i_row = rsiq_buf[te_type] + ((mo_i*BasisSet.num_ao + sq_fao) * nr * ns);
#if !USE_BLAS
			for(qrs=0;qrs<nq*nr*ns;qrs++,i_row++,rspq_ptr++) {
			  (*i_row) += temp * (*rspq_ptr);
			}
#else
			C_DAXPY(nq*nr*ns,temp,rspq_ptr,1,i_row,1);
			rspq_ptr += nq*nr*ns;
#endif
		      }
		    }
		  else
		    for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		      mo_vec = MOInfo.scf_evec_occ[0][mo_i+imin];
		      rspq_ptr = data[te_type];
		      for(p=0,p_abs=sp_fao;p<np;p++,p_abs++) {
			temp = mo_vec[p_abs];
			i_row = rsiq_buf[te_type] + ((mo_i*BasisSet.num_ao + sq_fao) * nr * ns);
#if !USE_BLAS
			for(qrs=0;qrs<nq*nr*ns;qrs++,i_row++,rspq_ptr++) {
			  (*i_row) += temp * (*rspq_ptr);
			}
#else
			C_DAXPY(nq*nr*ns,temp,rspq_ptr,1,i_row,1);
			rspq_ptr += nq*nr*ns;
#endif
		      }
		    }
		}
		timer_off("Step 1");
	      
	      } /* end of computing "petit" list - end of P,Q loop */
	    } /* end of "unique" P,Q loop */
	  
	  /*--- step 2 of the transfromation ---*/
	  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	    timer_on("Post1Swap");
	    iq_row = ijkl_to_klij(rsiq_buf[te_type],ibatch_length*BasisSet.num_ao,nr*ns);
	    timer_off("Post1Swap");
	    rsi_row = iq_row;
	    for(r=0;r<nr;r++) {
	      for(s=0;s<ns;s++) {
		  /*--- Can be done as a matrix multiply now ---*/
		  timer_on("Step 2");
#if !USE_BLAS
		  for(mo_i=0;mo_i<ibatch_length;mo_i++,rsi_row+=BasisSet.num_ao) {
		    for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
		      mo_vec = MOInfo.scf_evec[0][mo_x];
		      temp = ix_buf[mo_i][mo_x];
		      for(q_abs=0;q_abs<BasisSet.num_ao;q_abs++) {
			temp += mo_vec[q_abs] * rsi_row[q_abs];
		      }
		      ix_buf[mo_i][mo_x] = temp;
		    }
		  }
#else
		  C_DGEMM('n','t',ibatch_length,MOInfo.num_mo,BasisSet.num_ao,1.0,
			  rsi_row,BasisSet.num_ao,MOInfo.scf_evec[0][0],BasisSet.num_ao,
			  0.0,ix_buf[0],MOInfo.num_mo);
		  rsi_row += BasisSet.num_ao*ibatch_length;
#endif
		  timer_off("Step 2");

		  /*--- step 3 of the transformation ---*/
		  timer_on("Step 3");
		  r_abs = r + sr_fao;
		  s_abs = s + ss_fao;
		  for(mo_i=0;mo_i<ibatch_length;mo_i++) {
		    i_row = ix_buf[mo_i];
		    for(mo_j=0;mo_j<MOInfo.nactdocc;mo_j++) {
		      jsi_row = jsix_buf[te_type] + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
		      temp = MOInfo.scf_evec_occ[0][mo_j+jmin][r_abs];
		      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
			jsi_row[mo_x] += temp * i_row[mo_x];
		      }
		      if (usi != usj) {
			jsi_row = jsix_buf[te_type] + ((mo_j * BasisSet.num_ao + r_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
			temp = MOInfo.scf_evec_occ[0][mo_j+jmin][s_abs];
			if (te_type == 2)   /*--- [r12,T1] integral - negative sign ---*/
			  temp *= (-1.0);
			for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
			  jsi_row[mo_x] += temp * i_row[mo_x];
			}
		      }
		    }
		  }
		  memset(ix_buf[0],0,ibatch_length*MOInfo.num_mo*sizeof(double));
		  timer_off("Step 3");
	      }
	    }
	    memset(rsiq_buf[te_type],0,nr*ns*ibatch_length*BasisSet.num_ao*sizeof(double));
	  }
	  
	  } /* end of R,S loop */
      } /* end of "unique" R,S loop */

    /*--- step 4 of the transformation ---*/
    timer_on("Step 4");
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<MOInfo.nactdocc;mo_j++) {
	  for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	    jyi_row = jyix_buf[te_type] + ((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo;
	    for(s_abs=0;s_abs<BasisSet.num_ao;s_abs++) {
	      jsi_row = jsix_buf[te_type] + ((mo_j * BasisSet.num_ao + s_abs) * ibatch_length + mo_i) * MOInfo.num_mo;
	      temp = MOInfo.scf_evec[0][mo_y][s_abs];
	      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
		jyi_row[mo_x] += temp * jsi_row[mo_x];
	      }
	    }
	  }
	}
      }
    }
    timer_off("Step 4");

    /*--- Finish the symmetrization step - zero out non-totally symmetric integrals in Abelian case */
/*    if (Symmetry.nirreps > 1)
      for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	  for(mo_j=0;mo_j<MOInfo.nactdocc;mo_j++) {
	    for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	      for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
		if ((MOInfo.mo2symblk_occ[0][mo_i+imin] ^ MOInfo.mo2symblk_occ[0][mo_j+jmin]) ^
		    (MOInfo.mo2symblk[mo_x] ^ MOInfo.mo2symblk[mo_y]))
		  jyix_buf[te_type][((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo + mo_x] = 0.0;
	      }
	    }
	  }
	}
      }*/

    /*--- Print them out for now ---*/
/*    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
      fprintf(outfile,"  Transformed integrals of type %d\n",te_type);
      for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	for(mo_j=0;mo_j<MOInfo.nactdocc;mo_j++) {
	  for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
	    for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
	      temp = jyix_buf[te_type][((mo_j * MOInfo.num_mo + mo_y) * ibatch_length + mo_i) * MOInfo.num_mo + mo_x];
	      if (fabs(temp) > ZERO) {
		if (te_type < 2)
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_j+jmin,mo_y,mo_i+imin,mo_x,
			mo_j*MOInfo.num_mo+mo_y,
			mo_i*MOInfo.num_mo+mo_x,
			temp);
		else
		fprintf(outfile,"<%d %d %d %d [%d] [%d] = %20.10lf\n",
			mo_i+imin,mo_x,mo_j+jmin,mo_y,
			mo_i*MOInfo.num_mo+mo_x,
			mo_j*MOInfo.num_mo+mo_y,
			temp);
	      }
	    }
	  }
	}
      }
    }*/

    /*--------------------------------------------------------
      Dump all fully transformed integrals to disk including
      the ones corresponding to non-active orbitals. That way
      it's easy to recompute MP2-R12 energy with and without
      frozen core later on. Zero out
      non-symmetrical integrals (what's left of the Pitzer's
      equal contribution theorem in Abelian case).
     --------------------------------------------------------*/
    for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	/*-------------------------------------------
	  Swap jy and ix again, this time to satisfy
	  IWL routines
	 -------------------------------------------*/
	ixjy_buf[te_type] = ijkl_to_klij(jyix_buf[te_type],
				     MOInfo.nactdocc*MOInfo.num_mo,
				     ibatch_length*MOInfo.num_mo);
	/*--------------------------------------------------------------------
	  Here's the tricky part. IWL routine wants all integrals with common
	  first two indices (i and x) and the same symmetry (jsym) of the
	  third index (j) in a buffer (jy_buf) of dimensions jrange x yrange,
	  where jrange is the number of DOCCs in jsym block, and yrange is
	  the number MOs in ysym block. We have to do some serious reindexing
	  here since we transformed only the integrals with active DOCC i
	  and j, etc.
	 --------------------------------------------------------------------*/
	for(mo_i=0;mo_i<ibatch_length;mo_i++) {
	  i = act2fullQTS[mo_i+imin];       /*--- mo_i+imin is the index in QTS "active" indexing scheme
					          we need i to be in QTS "full" indexing scheme, that's
						  what IWL and MP2R12 code expect ---*/
	  isym = MOInfo.mo2symblk_occ[0][mo_i+imin];
	  for(mo_x=0;mo_x<MOInfo.num_mo;mo_x++) {
	    x = mo_x;    /*--- The second index is a Pitzer index, that's fine ---*/
	    xsym = MOInfo.mo2symblk[x];
	    /*--- Put all integrals with common i and x into a buffer ---*/
	    for(jsym = 0;jsym<Symmetry.nirreps;jsym++) {
	      /*--- if there're any active DOCCs in this jsym - proceed ---*/
	      if (MOInfo.clsdpi[jsym] - MOInfo.frozen_docc[jsym] != 0) {
		memset(jy_buf[0],0,MOInfo.ndocc*MOInfo.num_mo*sizeof(double));
		for(mo_j=0;mo_j<MOInfo.nactdocc;mo_j++) {
		  /*--- if this j is in jsym block - proceed ---*/
		  if (jsym == MOInfo.mo2symblk_occ[0][mo_j+jmin]) {
		    j = act2fullQTS[mo_j+jmin];    /*--- Again, get the "full" QTS index ---*/
		    for(mo_y=0;mo_y<MOInfo.num_mo;mo_y++) {
		      y = mo_y;                    /*--- Again, the Pitzer index here is what we need ---*/
		      ysym = MOInfo.mo2symblk[y];
		      /*--- Skip this indegral if it's non-totally symmetric -
			    Pitzer's contribution theorem in Abelian case ---*/
		      if ((isym ^ jsym) ^ (xsym ^ ysym))
			continue;
		      /*--- find the integral in ixjy_buf ---*/
		      m = (((mo_i*MOInfo.num_mo + mo_x)*MOInfo.nactdocc + mo_j)*MOInfo.num_mo + mo_y);
		      /*--- Now, the buffer is expected to be of jrange by yrange, hence we
			have to reindex to know where to put the integral. The required
			indices are simply relative indices in symmetry blocks - so far we've had
			absolute indices - hence have to subtract offsets
			first[jsym] - the offset for jsym in Pitzer scheme
			occ[fstocc[jsym]] - maps Pitzer offset index to "full" QTS DOCC index
			first[ysym] - the offset for ysym in Pitzer scheme ---*/
		      jy_buf[j-occ[first[jsym]]][y-first[ysym]] = ixjy_buf[te_type][m];
		    }
		  }
		}
		iwl_buf_wrt_mp2r12a(&MOBuf[te_type], i, x, ioff3[i]+x,
				    isym^xsym, jy_buf, jsym,
				    fstocc, lstocc, first, last,
				    occ,
				    (te_type == 2) ? 0 : 1,    /*--- No bra-ket permutational symmetry
								     for [r12,T2] integrals ---*/
				    ioff3, 0, outfile);
	      }
	    }
	  }
	}
    }

    if (ibatch < num_ibatch-1) {
      for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
	memset(jsix_buf[te_type],0,MOInfo.nactdocc*BasisSet.num_ao*ibatch_length*MOInfo.num_mo*sizeof(double));
	memset(jyix_buf[te_type],0,MOInfo.nactdocc*MOInfo.num_mo*ibatch_length*MOInfo.num_mo*sizeof(double));
      }
    }

  } /* End of "I"-loop */

  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    iwl_buf_flush(&MOBuf[te_type], 1); 
    iwl_buf_close(&MOBuf[te_type], 1);
  }
  fprintf(outfile,"  Transformation finished. Use the MP2R12 program to compute MP2-R12 energy.\n");
  fprintf(outfile,"  WARNING: Please, use the same FROZEN_DOCC vector with MP2-R12 as here,\n");
  fprintf(outfile,"           otherwise you will get meaningless results.\n\n");
  
  /*---------
    Clean-up
   ---------*/
  for(te_type=0;te_type<NUM_TE_TYPES-1;te_type++) {
    free(jyix_buf[te_type]);
    free(jsix_buf[te_type]);
    free(rsiq_buf[te_type]);
  }
  free_block(ix_buf);
  free_block(jy_buf);
  free_libr12(&Libr12);
  free_fjt_table(&fjt_table);
#ifdef NONDOUBLE_INTS
  for(i=0;i<NUM_TE_TYPES;i++)
    free(raw_data[i]);
#endif
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
  free_fjt();
  timer_done();

  return;
}




void make_transqt_arrays(int **first, int **last, int **fstocc, int **lstocc, int **occ, int **act2fullQTS, int **ioff3)
{
  int h, i, offset, count;
  int first_offset, last_offset;
  int p,q,row,col;

  /*
     Construct first and last index arrays: this defines the first
     absolute orbital index (Pitzer ordering) and last absolute orbital
     index for each irrep.  When there are no orbitals for an irrep, the
     value is -1 for first[] and -2 for last[].  Note that there must be
     orbitals in the first irrep (i.e. totally symmetric) for this to work.
  */
  *first = init_int_array(Symmetry.nirreps);
  *last = init_int_array(Symmetry.nirreps);
  for(h=0; h < Symmetry.nirreps; h++) {
      (*first)[h] = -1;
      (*last)[h] = -2;
  }
  first_offset = 0;
  last_offset = MOInfo.orbspi[0] - 1; 
  (*first)[0] = first_offset;
  (*last)[0] = last_offset;
  for(h=1; h < Symmetry.nirreps; h++) {
      first_offset += MOInfo.orbspi[h-1];
      last_offset += MOInfo.orbspi[h];
      if(MOInfo.orbspi[h]) {
          (*first)[h] = first_offset;
          (*last)[h] = last_offset;
        }
    }
  
  /* fstocc[] and lstocc[] supply the first and last orbital indices (Pitzer
     ordering) for the occupied orbitals in each irrep. */
  *fstocc = init_int_array(Symmetry.nirreps);
  *lstocc = init_int_array(Symmetry.nirreps);
  for(h=0; h < Symmetry.nirreps; h++) {
      (*fstocc)[h] = -1;
      (*lstocc)[h] = -2;
    }
  first_offset = 0;
  last_offset = MOInfo.clsdpi[0] - 1;
  (*fstocc)[0] = first_offset;
  (*lstocc)[0] = last_offset;
  for(h=1; h < Symmetry.nirreps; h++) {
      first_offset += MOInfo.orbspi[h-1];
      last_offset += MOInfo.virtpi[h-1]+MOInfo.orbspi[h]-MOInfo.virtpi[h];
      if(MOInfo.clsdpi[h]) {
          (*fstocc)[h] = first_offset;
          (*lstocc)[h] = last_offset;
        }
    }

  /* Construct occupied Pitzer -> QTS ordering arrays for
   occupied (occ[]) orbitals */
  *occ = init_int_array(MOInfo.num_mo);
  for(i=0; i< MOInfo.num_mo; i++) {
      (*occ)[i] = -1;
    }
  
  offset = 0;
  count=0;
  for(h=0; h < Symmetry.nirreps; h++) {
      if(h)
          offset += MOInfo.orbspi[h-1];
      for(i=offset; i < (offset+MOInfo.clsdpi[h]); i++) {
          (*occ)[i] = count++;
        }
    }

  /* Construct active -> full QTS ordering array for occupied orbitals */
  *act2fullQTS = init_int_array(MOInfo.ndocc);
  offset = 0;
  count=0;
  for(h=0; h < Symmetry.nirreps; h++) {
    for(i=0; i < MOInfo.frozen_docc[h]; i++) {
      (*act2fullQTS)[count] = offset+i;
      count++;
    }
    offset += MOInfo.clsdpi[h];
  }
  offset = 0;
  for(h=0; h < Symmetry.nirreps; h++) {
    for(i=MOInfo.frozen_docc[h]; i < MOInfo.clsdpi[h]; i++) {
      (*act2fullQTS)[count] = offset+i;
      count++;
    }
    offset += MOInfo.clsdpi[h];
  }
  
  /* Generate ioff3 array.  This array gives the row offset for an
     ndocc x nmo matrix */
  *ioff3 = init_int_array(MOInfo.ndocc);
  for(i=0; i < MOInfo.ndocc; i++) {
      (*ioff3)[i] = i*MOInfo.num_mo;
    }

  return;
}
