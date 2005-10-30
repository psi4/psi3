#include <vector>

#include<math.h>
#include<stdio.h>
#include<string.h>
#include<memory.h>
#include<stdlib.h>
#include<libipv1/ip_lib.h>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>

#include <libint/libint.h>
#include<libderiv/libderiv.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include"schwartz.h"
#include"deriv1_quartet_data.h"
#include"iwl_tebuf.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
  #include"taylor_fm_eval.h"
#else
  #include"int_fjt.h"
#endif

void te_deriv1_ints()
{
  const double toler = UserOptions.cutoff;
  const int num_coords = 3*Molecule.num_atoms;

  /*--- ASCII file to print integrals ---*/
  FILE *d1eriout;

  /*--- Various data structures ---*/
  struct iwlbuf* D1ERIOUT;               /* IWL buffer for target integrals */
  struct tebuf** tot_data;                /* accum. for contracted integrals */
  struct shell_pair *sp_ij, *sp_kl;
  struct unique_shell_pair *usp_ij,*usp_kl;

  Libderiv_t Libderiv;
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;
#endif

  int total_te_count = 0;
  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  int i, j, k, l, m, ii, jj, kk, ll;
  int si, sj, sk, sl ;
  int sii, sjj, skk, sll , slll;
  int pi, pj, pk, pl;
  int max_pj, max_pl;
  int upk, num_unique_pk;
  int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
  int *sj_arr, *sk_arr, *sl_arr;
  int *sj_fbf_arr, *sk_fbf_arr, *sl_fbf_arr;
  int usi,usj,usk,usl;
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

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;

  double so_int;
  double AB2, CD2;
  double *data;                 /* pointer to the transformed normalized target quartet of integrals */
  double *puream_data;
  double **plist_data;
  double pkblock_end_value = 0.0;
  double temp;

  /*---------------
    Initialization
   ---------------*/
#if PRINT
  eriout = fopen("d1eriout.dat","w");
#endif
  D1ERIOUT = new struct iwlbuf[num_coords];
  for(int c=0; c<num_coords; c++)
    iwl_buf_init(&D1ERIOUT[c], IOUnits.itapD1ERI_SO+c, toler, 0, 0);
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+DERIV_LVL,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am*4+DERIV_LVL);
  init_fjt_table(&fjt_table);
#endif
  init_libderiv_base();
  
  /*-------------------------
    Allocate data structures
   -------------------------*/
  max_cart_class_size = (ioff[BasisSet.max_am])*
                        (ioff[BasisSet.max_am])*
                        (ioff[BasisSet.max_am])*
                        (ioff[BasisSet.max_am]);
  max_num_unique_quartets = Symmetry.max_stab_index*
                            Symmetry.max_stab_index*
                            Symmetry.max_stab_index;
  
  // Final integrals are stored here
  tot_data = new struct tebuf*[num_coords];
  for(int c=0; c<num_coords; c++) {
    tot_data[c] = new struct tebuf[max_num_unique_quartets*max_cart_class_size];
    memset(tot_data[c], 0, (max_num_unique_quartets*max_cart_class_size)*sizeof(struct tebuf));
  }
  // These arrays are used to hold cartesian and spherical harmonics integrals
  double* cart_ints[12];
  for(int i=0; i<12; i++)
    cart_ints[i] = new double[max_cart_class_size];
  double* sph_ints[12];
  if (BasisSet.puream) {
    for(int i=0; i<12; i++)
      sph_ints[i] = new double[max_cart_class_size];
  }
  else {
    for(int i=0; i<12; i++)
      sph_ints[i] = cart_ints[i];
  }
  
  sj_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sk_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  sl_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  if (Symmetry.nirreps > 1) {
    if (BasisSet.puream)
      max_class_size = (2*BasisSet.max_am-1)*(2*BasisSet.max_am-1)*(2*BasisSet.max_am-1)*(2*BasisSet.max_am-1);
    else
      max_class_size = max_cart_class_size;
    plist_data = block_matrix(max_num_unique_quartets,max_class_size);
    sj_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sk_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
    sl_fbf_arr = (int *)malloc(sizeof(int)*max_num_unique_quartets);
  }
  
  max_num_prim_comb = (BasisSet.max_num_prims*
                       BasisSet.max_num_prims)*
                      (BasisSet.max_num_prims*
                       BasisSet.max_num_prims);
  init_libderiv1(&Libderiv, BasisSet.max_am-1, max_num_prim_comb, max_class_size);

  /*------------------------------------
    generate all unique shell quartets
   ------------------------------------*/
  if(UserOptions.print_lvl >= PRINT_DEBUG) {
    fprintf(outfile,"  -Electron repulsion integrals:\n\n");
  }
  
  for (usi=0; usi<Symmetry.num_unique_shells; usi++)
    for (usj=0; usj<=usi; usj++)
      for (usk=0; usk<Symmetry.num_unique_shells; usk++) {
        int usl_max = (usi==usk ? usj : usk);
	for (usl=0; usl<=usl_max; usl++){

          /* place in "ascending" angular mom-
          my simple way of optimizing PHG recursion (VRR) */
          /* these first two are good for the HRR */
          if(BasisSet.shells[Symmetry.us2s[usi]].am < BasisSet.shells[Symmetry.us2s[usj]].am){
            dum = usi;
            usi = usj;
            usj = dum;
          }
          if(BasisSet.shells[Symmetry.us2s[usk]].am < BasisSet.shells[Symmetry.us2s[usl]].am){
            dum = usk;
            usk = usl;
            usl = dum;
          }
          /* this should be /good/ for the VRR */
          if(BasisSet.shells[Symmetry.us2s[usi]].am + BasisSet.shells[Symmetry.us2s[usj]].am >
          BasisSet.shells[Symmetry.us2s[usk]].am + BasisSet.shells[Symmetry.us2s[usl]].am){
            dum = usi;
            usi = usk;
            usk = dum;
            dum = usj;
            usj = usl;
            usl = dum;
          }
          
          si = Symmetry.us2s[usi];
          sjj = Symmetry.us2s[usj];
          skk = Symmetry.us2s[usk];
          sll = Symmetry.us2s[usl];
          
          if (Symmetry.nirreps > 1) { /*--- Non-C1 symmetry case ---*/
            /*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
            usp_ij = &(Symmetry.us_pairs[usi][usj]);
            usp_kl = &(Symmetry.us_pairs[usk][usl]);
            stab_i = Symmetry.atom_positions[BasisSet.shells[si].center-1];
            stab_j = Symmetry.atom_positions[BasisSet.shells[sjj].center-1];
            stab_k = Symmetry.atom_positions[BasisSet.shells[skk].center-1];
            stab_l = Symmetry.atom_positions[BasisSet.shells[sll].center-1];
            stab_ij = Symmetry.GnG[stab_i][stab_j];
            stab_kl = Symmetry.GnG[stab_k][stab_l];
            R_list = Symmetry.dcr[stab_i][stab_j];
            S_list = Symmetry.dcr[stab_k][stab_l];
            T_list = Symmetry.dcr[stab_ij][stab_kl];
            lambda_T = Symmetry.nirreps/Symmetry.dcr_deg[stab_ij][stab_kl];
            ni = (BasisSet.puream ? 2*BasisSet.shells[si].am - 1 : ioff[BasisSet.shells[si].am]);
            nj = (BasisSet.puream ? 2*BasisSet.shells[sjj].am - 1 : ioff[BasisSet.shells[sjj].am]);
            nk = (BasisSet.puream ? 2*BasisSet.shells[skk].am - 1 : ioff[BasisSet.shells[skk].am]);
            nl = (BasisSet.puream ? 2*BasisSet.shells[sll].am - 1 : ioff[BasisSet.shells[sll].am]);
            class_size = ni*nj*nk*nl;
            
            memset(sj_arr,0,sizeof(int)*max_num_unique_quartets);
            memset(sk_arr,0,sizeof(int)*max_num_unique_quartets);
            memset(sl_arr,0,sizeof(int)*max_num_unique_quartets);
            memset(sj_fbf_arr,0,sizeof(int)*max_num_unique_quartets);
            memset(sk_fbf_arr,0,sizeof(int)*max_num_unique_quartets);
            memset(sl_fbf_arr,0,sizeof(int)*max_num_unique_quartets);
            count = 0;
            for(dcr_ij=0;dcr_ij<Symmetry.dcr_dim[stab_i][stab_j];dcr_ij++){
              R = R_list[dcr_ij];
              sj = BasisSet.shells[sjj].trans_vec[R]-1;
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
                    sj_fbf_arr[count] = BasisSet.shells[sj].fbf-1;
                    sk_fbf_arr[count] = BasisSet.shells[sk].fbf-1;
                    sl_fbf_arr[count] = BasisSet.shells[sl].fbf-1;
                    count++;
                  }
                }
              }
            } /* petite list is ready to be used */
            num_unique_quartets = count;
          }
          else { /*--- C1 symmetry case ---*/
            num_unique_quartets = 1;
            sj_arr[0] = usj;
            sk_arr[0] = usk;
            sl_arr[0] = usl;
            ni = (BasisSet.puream ? 2*BasisSet.shells[usi].am - 1 : ioff[BasisSet.shells[usi].am]);
            nj = (BasisSet.puream ? 2*BasisSet.shells[usj].am - 1 : ioff[BasisSet.shells[usj].am]);
            nk = (BasisSet.puream ? 2*BasisSet.shells[usk].am - 1 : ioff[BasisSet.shells[usk].am]);
            nl = (BasisSet.puream ? 2*BasisSet.shells[usl].am - 1 : ioff[BasisSet.shells[usl].am]);
            ioffset = BasisSet.shells[usi].fbf - 1;
            joffset = BasisSet.shells[usj].fbf - 1;
            koffset = BasisSet.shells[usk].fbf - 1;
            loffset = BasisSet.shells[usl].fbf - 1;
          }
          
          np_i = BasisSet.shells[si].n_prims;
          np_j = BasisSet.shells[sjj].n_prims;
          np_k = BasisSet.shells[skk].n_prims;
          np_l = BasisSet.shells[sll].n_prims;
          
          orig_am[0] = BasisSet.shells[si].am-1;
          orig_am[1] = BasisSet.shells[sjj].am-1;
          orig_am[2] = BasisSet.shells[skk].am-1;
          orig_am[3] = BasisSet.shells[sll].am-1;
          am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];
          

          /*----------------------------------
          Compute the nonredundant quartets
          ----------------------------------*/
          for(plquartet=0;plquartet<num_unique_quartets;plquartet++) {
            sj = sj_arr[plquartet];
            sk = sk_arr[plquartet];
            sl = sl_arr[plquartet];
            
            sp_ij = &(BasisSet.shell_pairs[si][sj]);
            sp_kl = &(BasisSet.shell_pairs[sk][sl]);
            
            Libderiv.AB[0] = sp_ij->AB[0];
            Libderiv.AB[1] = sp_ij->AB[1];
            Libderiv.AB[2] = sp_ij->AB[2];
            Libderiv.CD[0] = sp_kl->AB[0];
            Libderiv.CD[1] = sp_kl->AB[1];
            Libderiv.CD[2] = sp_kl->AB[2];
            
            AB2 = Libderiv.AB[0]*Libderiv.AB[0]+Libderiv.AB[1]*Libderiv.AB[1]+Libderiv.AB[2]*Libderiv.AB[2];
            CD2 = Libderiv.CD[0]*Libderiv.CD[0]+Libderiv.CD[1]*Libderiv.CD[1]+Libderiv.CD[2]*Libderiv.CD[2];
            
	    /*--- Compute data for primitive quartets here ---*/
	    num_prim_comb = 0;
            const double pfac = lambda_T;
	    for (pi = 0; pi < np_i; pi++) {
	      max_pj = (si == sj) ? pi+1 : np_j;
	      for (pj = 0; pj < max_pj; pj++) {
		m = (1 + (si == sj && pi != pj));
		for (pk = 0; pk < np_k; pk++) {
		  max_pl = (sk == sl) ? pk+1 : np_l;
		  for (pl = 0; pl < max_pl; pl++){
		    n = m * (1 + (sk == sl && pk != pl));
#ifdef USE_TAYLOR_FM
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					NULL, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
#else
		    deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
					&fjt_table, AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
#endif		    
		  }
		}
	      }
	    }
            
            /*--- Compute the derivative integrals ---*/
            size = ioff[BasisSet.shells[si].am]*ioff[BasisSet.shells[sj].am]*
                   ioff[BasisSet.shells[sk].am]*ioff[BasisSet.shells[sl].am];
            build_deriv1_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libderiv, num_prim_comb);
            
            for(int c=0; c<12; c++) {
              if (c < 3 || c > 5) {
#ifdef NONDOUBLE_INTS
                for(int j=0;j<size;j++)
                  cart_ints[c][j] = (double) Libderiv.ABCD[c][j];
#else
                cart_ints[c] = Libderiv.ABCD[c];
#endif
              }
            }
            // reconstruct integrals using translational invariance condition
            for(int j=0;j<size;j++) {
              cart_ints[3][j] = -(cart_ints[0][j] + cart_ints[6][j] + cart_ints[9][j]);
              cart_ints[4][j] = -(cart_ints[1][j] + cart_ints[7][j] + cart_ints[10][j]);
              cart_ints[5][j] = -(cart_ints[2][j] + cart_ints[8][j] + cart_ints[11][j]);
            }
            
            // determine all unique centers (0..3)
            std::vector<int> unique_center;
            std::vector<int> nuc;
            nuc.push_back(BasisSet.shells[si].center-1);
            nuc.push_back(BasisSet.shells[sj].center-1);
            nuc.push_back(BasisSet.shells[sk].center-1);
            nuc.push_back(BasisSet.shells[sl].center-1);
            unique_center.push_back(0);
            if (nuc[1] != nuc[0]) {
              unique_center.push_back(1);
            }
            if (nuc[2] != nuc[0] && nuc[2] != nuc[1]) {
              unique_center.push_back(2);
            }
            if (nuc[3] != nuc[0] && nuc[3] != nuc[1] && nuc[3] != nuc[2]) {
              unique_center.push_back(3);
            }
            const int num_unique_centers = unique_center.size();

            // compute total derivatives from partial derivatives
            for(int i=0; i<num_unique_centers; i++) {
              const int c = unique_center[i];
              for(int j=c+1; j<4; j++) {
                if (nuc[c] == nuc[j]) {
                  for(int xyz=0; xyz<3; xyz++) {
                    double* tot_der = cart_ints[3*c+xyz];
                    const double* part_der = cart_ints[3*j+xyz];
                    for(int k=0; k<size; k++)
                      tot_der[k] += part_der[k];
                  }
                }
              }
            }
            const int num_tot_der = num_unique_centers*3;
            
            for(int d=0; d<num_tot_der; d++) {
              const double* data = norm_quartet(cart_ints[d], puream_data, orig_am, BasisSet.puream);
              if (data == sph_ints[d])
                for(int k=0; k<size; k++)
                  sph_ints[d][k] = data[k];
            }
            /*--- copy data to plist_data to be used in the symmetrization step ---*/
            if (Symmetry.nirreps > 1)
              memcpy(plist_data[plquartet],data,sizeof(double)*class_size);
            
          } /* end of computing "petit" list */

	    num = 0;
	    if (Symmetry.nirreps > 1) { /*--- Non-C1 case ---*/
	      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
		  /*---
		    npi_ij - number of pairs of SOs arising from the ij pair of unique shells
		             whose direct product transforms as irrep
		    check to see if for given irrep npi_ij and npi_kl are non-zero, i.e.
		    if there's a combination of SO's that thansforms as totally symm. irrep
		    ---*/
		if ((npi_ij = usp_ij->SOpair_npi[irrep]) && (npi_kl = usp_kl->SOpair_npi[irrep]))
		  for(ij=0;ij<npi_ij;ij++) { /*--- Loop over SO pairs from usij ---*/
		    i = usp_ij->SOpair_bf_i[irrep][ij];       /*--- Basis function type to which this SO from usi corresponds,
								    e.g. px, dxz, etc. ---*/
		    j = usp_ij->SOpair_bf_j[irrep][ij];
		    so_i = usp_ij->SOpair_so_i[irrep][ij];    /*--- Absolute index of this SO from usi ---*/
		    so_j = usp_ij->SOpair_so_j[irrep][ij];
		    ind_offset = (i*nj + j)*nk*nl;
		    for(kl=0;kl<npi_kl;kl++) { /*--- Loop over SO pairs from uskl ---*/
		      k = usp_kl->SOpair_bf_i[irrep][kl];
		      l = usp_kl->SOpair_bf_j[irrep][kl];
		      so_k = usp_kl->SOpair_so_i[irrep][kl];
		      so_l = usp_kl->SOpair_so_j[irrep][kl];
                      
		      index = ind_offset + k*nl + l;          /*--- Based on basis function types figure
							 	    out which particular integral from these
								    quartets contributes ---*/
		      so_int = 0.0;
		      for(s=0;s<num_unique_quartets;s++){  /*--- Sum over petite list quartets ---*/
			ioffset = bf_i + i;      /*--- Absolute basis function index ---*/
			bf_j = sj_fbf_arr[s]+j;
			bf_k = sk_fbf_arr[s]+k;
			bf_l = sl_fbf_arr[s]+l;
			so_int += Symmetry.usotao[so_i][ioffset]*
				  Symmetry.usotao[so_j][bf_j]*
				  Symmetry.usotao[so_k][bf_k]*
				  Symmetry.usotao[so_l][bf_l]*
				  plist_data[s][index];
		      }
                    }
                  }
	      }
	    }
	    else { /*--- C1 symmetry ---*/
              iimax = ni - 1;
              kkmax = nk - 1;
              for(ii=0; ii <= iimax; ii++){
                jjmax = (usi==usj) ? ii : nj - 1;
                for(jj=0; jj <= jjmax; jj++){
                  for(kk=0; kk <= kkmax; kk++){
                    llmax = (usk==usl) ? kk : nl - 1;
                    for(ll=0; ll <= llmax; ll++){
                      index = ll+nl*(kk+nk*(jj+nj*ii));
                      if(fabs(data[index])>toler){
                        i = ii + ioffset;
                        j = jj + joffset;
                        k = kk + koffset;
                        l = ll + loffset;
                        if (i < j)
                          SWAP(i,j);
                        if (k < l)
                          SWAP(k,l);
                        if ((i < k) || (i == k && j < l)) {
                          SWAP(i,k);
                          SWAP(j,l);
                        }
                        /*tot_data[num].i = (short int) i;
                        tot_data[num].j = (short int) j;
                        tot_data[num].k = (short int) k;
                        tot_data[num].l = (short int) l;
                        tot_data[num].val = data[index];*/
                        num++;
                      }
                    }
                  }
                }
              }
            }

	    if (num) { /* Let's see if we need to write out something */
	      total_te_count += num;
              // iwl_buf_wrt_struct_nocut(&ED1RIOUT[c], tot_data[c], num);
	    }

	    if(UserOptions.print_lvl >= PRINT_DEBUG) {
	      /* fprintf(outfile,"  -Electron repulsion integrals:\n\n"); */
	    }  
	    
      }} /* end getting unique shell combination */
  
      for(int c=0; c<num_coords; c++) {
        iwl_buf_flush(&D1ERIOUT[c], 1);
        iwl_buf_close(&D1ERIOUT[c], 1);
      }

      
      fprintf(outfile,"\n    Wrote %d two-electron integrals to IWL file %2d\n\n",total_te_count,IOUnits.itap33);
      
      /*---------
      Clean-up
      ---------*/
      free_libderiv(&Libderiv);
      free(sj_arr);
      free(sk_arr);
      free(sl_arr);
      if (Symmetry.nirreps > 1) {
        free_block(plist_data);
        free(sj_fbf_arr);
        free(sk_fbf_arr);
        free(sl_fbf_arr);
      }
      #ifdef USE_TAYLOR_FM
      free_Taylor_Fm_Eval();
      #else
      free_fjt_table(&fjt_table);
      free_fjt();
      #endif

  return;
}

