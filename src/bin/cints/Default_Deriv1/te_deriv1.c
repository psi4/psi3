#include <math.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <iwl.h>
#include <libciomr.h>
#include <libint.h>
#include <libderiv.h>
#include "defines.h"
#define EXTERN
#include "global.h"
#include "int_fjt.h"
#include "deriv1_quartet_data.h"
#include "small_fns.h"


void te_deriv1()
{

  /*--- Various data structures ---*/
  struct iwlbuf TPDM;                 /* IWL buffer for two-pdm matrix elements */
  struct shell_pair *sp_ij, *sp_kl;

  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count ;
  int read_dens = strcmp(UserOptions.wfn,"SCF");
  int dum;
  int n, num;
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  register int i, j, k, l, m, ii, jj, kk, ll;
  register int si, sj, sk, sl ;
  register int sii, sjj, skk, sll, slll;
  register int pi, pj, pk, pl ;
  int max_pj, max_pl;
  register int pii, pjj, pkk, pll ;
  int switch_ij, switch_kl, switch_ijkl;
  int center_i, center_j, center_k, center_l;

  int class_size;
  int max_class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl, quartet_size;

  int si_fao, sj_fao, sk_fao, sl_fao;
  int sii_fao, sjj_fao, skk_fao, sll_fao;
  int ao_i, imax, ao_j, jmax, ao_k, kmax, ao_l, lmax;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p, max_num_prim_comb;

  int buf_offset, buf_4offset, buf_size, last_buf;
  int quartet_done, offset;

  int mosh_i, mosh_j;

  double AB2, CD2;
  double *FourInd;
  double **grad_te;
  double pfac;
  double temp;
  double alpha, beta;
  double **dens_i, **dens_j;
  double ddax, dday, ddaz, ddbx, ddby, ddbz,
         ddcx, ddcy, ddcz, dddx, dddy, dddz;


  /*---------------
    Initialization
   ---------------*/
  if (strcmp(UserOptions.wfn,"SCF")) {
    iwl_buf_init(&TPDM, IOUnits.itapG, 0.0, 1, 1);
    buf_offset = 0;
    buf_4offset = 0;
    buf_size = TPDM.inbuf;
  }
  int_initialize_fjt(BasisSet.max_am*4+DERIV_LVL);
  init_libderiv();

  
  /*-------------------------
    Allocate data structures
   -------------------------*/
  max_cart_class_size = ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am]*ioff[BasisSet.max_am];
  max_class_size = max_cart_class_size;

  int_stack = (double *) malloc(STACK_SIZE*sizeof(double));
  zero_stack = init_array(max_cart_class_size);
  max_num_prim_comb = (BasisSet.max_num_prims*BasisSet.max_num_prims)*
		      (BasisSet.max_num_prims*BasisSet.max_num_prims);
  Shell_Data = (prim_data *) malloc(max_num_prim_comb*sizeof(prim_data));
  FourInd = init_array(max_cart_class_size);

  grad_te = block_matrix(Molecule.num_atoms,3);

/*-------------------------------------------------
  generate all unique shell quartets with ordering
  suitable for building the PK-matrix
 -------------------------------------------------*/
  for (sii=0; sii<BasisSet.num_shells; sii++)
    for (sjj=0; sjj<=sii; sjj++)
      for (skk=0; skk<=sii; skk++)
	for (sll=0; sll<= ((sii == skk) ? sjj : skk); sll++){

	    si = sii; sj = sjj; sk = skk; sl = sll;

	    /*--- Skip this quartet if all four centers are the same ---*/
	    if (BasisSet.shells[si].center == BasisSet.shells[sj].center &&
		BasisSet.shells[si].center == BasisSet.shells[sk].center &&
		BasisSet.shells[si].center == BasisSet.shells[sl].center &&
		DERIV_LVL == 1) {
              /*--- If reading in density - need to skip the appropriate block of that too ---*/
/*	      if (read_dens) {
		last_buf = TPDM.lastbuf;
		quartet_done = 0;
		do {
		  if (buf_offset < buf_size) {
		    i = TPDM.labels[buf_4offset]   - si_fao;
		    if (i < 0)
		      quartet_done = 1;
		    buf_offset++;
		    buf_4offset += 4;
		  }
		  else if (!last_buf) {
		    iwl_buf_fetch(&TPDM);
		    buf_offset = 0;
		    buf_4offset = 0;
		    last_buf = TPDM.lastbuf;
		  }
		  else {
		    punt(fpo,"  The last TPDM quartet not marked\n");
		  }
		} while (!quartet_done);
	      }*/

	      continue;
	    }

	    switch_ij = 0;
	    switch_kl = 0;
	    switch_ijkl = 0;
	    /* place in "ascending" angular mom-
	       my simple way of optimizing PHG recursion (VRR) */
	    /* these first two are good for the HRR */
	    if(BasisSet.shells[si].am < BasisSet.shells[sj].am){
	      dum = si;
	      si = sj;
	      sj = dum;
	      switch_ij = 1;
	    }
	    if(BasisSet.shells[sk].am < BasisSet.shells[sl].am){
	      dum = sk;
	      sk = sl;
	      sl = dum;
	      switch_kl = 1;
	    }
	    /* this should be /good/ for the VRR */
	    if(BasisSet.shells[si].am + BasisSet.shells[sj].am > BasisSet.shells[sk].am + BasisSet.shells[sl].am){
	      dum = si;
	      si = sk;
	      sk = dum;
	      dum = sj;
	      sj = sl;
	      sl = dum;
	      switch_ijkl = 1;
	    }

	    ni = ioff[BasisSet.shells[si].am];
	    nj = ioff[BasisSet.shells[sj].am];
	    nk = ioff[BasisSet.shells[sk].am];
	    nl = ioff[BasisSet.shells[sl].am];
	    quartet_size = ni*nj*nk*nl;

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
	    
	    AB[0] = sp_ij->AB[0];
	    AB[1] = sp_ij->AB[1];
	    AB[2] = sp_ij->AB[2];
	    CD[0] = sp_kl->AB[0];
	    CD[1] = sp_kl->AB[1];
	    CD[2] = sp_kl->AB[2];
	    
	    AB2 = AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2];
	    CD2 = CD[0]*CD[0]+CD[1]*CD[1]+CD[2]*CD[2];

	    /*-------------------------
	      Figure out the prefactor
	     -------------------------*/
	    pfac = 1.0;
	    if (si == sj)
	      pfac *= 0.5;
	    if (sk == sl)
	      pfac *= 0.5;
	    if (si == sk && sj == sl || si == sl && sj == sk)
	      pfac *= 0.5;
	    if (read_dens) { /*--- The factor of 8 needed for correlated densities ---*/
	      pfac *= 8.0;
	    }

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
		    deriv1_quartet_data(&(Shell_Data[num_prim_comb++]), AB2, CD2,
					sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
		    
		  }
		}
	      }
	    }

	    /*-------------
	      Form FourInd
	     -------------*/
	    if (read_dens) { /*--- Read in a shell quartet from disk ---*/
	      memset(FourInd,0,sizeof(double)*quartet_size);
	      last_buf = TPDM.lastbuf;
	      quartet_done = 0;
	      sii_fao = BasisSet.shells[sii].fao-1;
	      sjj_fao = BasisSet.shells[sjj].fao-1;
	      skk_fao = BasisSet.shells[skk].fao-1;
	      sll_fao = BasisSet.shells[sll].fao-1;
	      do {
		if (buf_offset < buf_size) {
		  i = TPDM.labels[buf_4offset]   - sii_fao;
		  if (i >= 0) {
		    j = TPDM.labels[buf_4offset+1] - sjj_fao;
		    k = TPDM.labels[buf_4offset+2] - skk_fao;
		    l = TPDM.labels[buf_4offset+3] - sll_fao;
		    if (switch_ij) {
		      dum = i;
		      i = j;
		      j = dum;
		    }
		    if (switch_kl) {
		      dum = k;
		      k = l;
		      l = dum;
		    }
		    if (switch_ijkl) {
		      dum = i;
		      i = k;
		      k = dum;
		      dum = j;
		      j = l;
		      l = dum;
		    }
		    offset = ((i*nj+j)*nk+k)*nl+l;
		    FourInd[offset] = TPDM.values[buf_offset]*
				      GTOs.bf_norm[orig_am[0]][i]*
				      GTOs.bf_norm[orig_am[1]][j]*
				      GTOs.bf_norm[orig_am[2]][k]*
				      GTOs.bf_norm[orig_am[3]][l];  
		  }
		  else
		    quartet_done = 1;
		  buf_offset++;
		  buf_4offset += 4;
		}
		else if (!last_buf) {
		  iwl_buf_fetch(&TPDM);
		  buf_offset = 0;
		  buf_4offset = 0;
		  last_buf = TPDM.lastbuf;
		}
		else {
		  punt("The last TPDM quartet not marked");
		}
	      } while (!quartet_done);
	    }
	    else { /*--- Compute a shell quartet of TPDM ---*/
	      si_fao = BasisSet.shells[si].fao-1;
	      sj_fao = BasisSet.shells[sj].fao-1;
	      sk_fao = BasisSet.shells[sk].fao-1;
	      sl_fao = BasisSet.shells[sl].fao-1;
	      if (UserOptions.reftype == rhf || UserOptions.reftype == uhf) {           /*--- RHF or UHF ---*/
		count = 0;
		for (ao_i = si_fao; ao_i < si_fao+ni; ao_i++)
		  for (ao_j = sj_fao; ao_j < sj_fao+nj; ao_j++)
		    for (ao_k = sk_fao; ao_k < sk_fao+nk; ao_k++)
		      for (ao_l = sl_fao; ao_l < sl_fao+nl; ao_l++) {
			FourInd[count] = (4.0*Dens[ao_i][ao_j]*Dens[ao_k][ao_l] -
					  Dens[ao_i][ao_k]*Dens[ao_j][ao_l] -
					  Dens[ao_i][ao_l]*Dens[ao_k][ao_j])*
					 GTOs.bf_norm[orig_am[0]][ao_i-si_fao]*
					 GTOs.bf_norm[orig_am[1]][ao_j-sj_fao]*
					 GTOs.bf_norm[orig_am[2]][ao_k-sk_fao]*
					 GTOs.bf_norm[orig_am[3]][ao_l-sl_fao];
			count++;
		      }
	      }
	      else {                     /*--- ROHF or TCSCF ---*/
		if (am)
		  memset((char *) FourInd, 0, sizeof(double)*quartet_size);
		else
		  FourInd[0] = 0.0;
		for(mosh_i=0;mosh_i<MOInfo.num_moshells;mosh_i++)
		  for(mosh_j=0;mosh_j<MOInfo.num_moshells;mosh_j++) {
		    count = 0;
		    dens_i = ShDens[mosh_i];
		    dens_j = ShDens[mosh_j];
		    alpha = 8.0*MOInfo.Alpha[mosh_i][mosh_j];
		    beta = 4.0*MOInfo.Beta[mosh_i][mosh_j];
		    for (ao_i = si_fao; ao_i < si_fao+ni; ao_i++)
		      for (ao_j = sj_fao; ao_j < sj_fao+nj; ao_j++)
			for (ao_k = sk_fao; ao_k < sk_fao+nk; ao_k++)
			  for (ao_l = sl_fao; ao_l < sl_fao+nl; ao_l++) {
			      FourInd[count] += (alpha*dens_i[ao_i][ao_j]*dens_j[ao_k][ao_l] +
						beta*(dens_i[ao_i][ao_k]*dens_j[ao_j][ao_l] +
						dens_i[ao_i][ao_l]*dens_j[ao_k][ao_j]));
			      count++;
			  }
		  }
		/*--- Normalize it ---*/
		count = 0;
		for (ao_i = 0; ao_i < ni; ao_i++)
		  for (ao_j = 0; ao_j < nj; ao_j++)
		    for (ao_k = 0; ao_k < nk; ao_k++)
		      for (ao_l = 0; ao_l < nl; ao_l++) {
			FourInd[count] *= GTOs.bf_norm[orig_am[0]][ao_i]*
					 GTOs.bf_norm[orig_am[1]][ao_j]*
					 GTOs.bf_norm[orig_am[2]][ao_k]*
					 GTOs.bf_norm[orig_am[3]][ao_l];
			count++;
		      }
	      }
	    }
	      
	    dbuild_abcd[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](num_prim_comb);

	    center_i = BasisSet.shells[si].center-1;
	    center_j = BasisSet.shells[sj].center-1;
	    center_k = BasisSet.shells[sk].center-1;
	    center_l = BasisSet.shells[sl].center-1;
	    
	    ddax = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddax += ABCD[0][k]*FourInd[k];
	    grad_te[center_i][0] += ddax;

	    dday = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dday += ABCD[1][k]*FourInd[k];
	    grad_te[center_i][1] += dday;

	    ddaz = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddaz += ABCD[2][k]*FourInd[k];
	    grad_te[center_i][2] += ddaz;

	    ddbx = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddbx += ABCD[3][k]*FourInd[k];
	    grad_te[center_j][0] += ddbx;

	    ddby = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddby += ABCD[4][k]*FourInd[k];
	    grad_te[center_j][1] += ddby;

	    ddbz = 0.0;
	    for(k=0;k<quartet_size;k++)
	      ddbz += ABCD[5][k]*FourInd[k];
	    grad_te[center_j][2] += ddbz;

	    dddx = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dddx += ABCD[9][k]*FourInd[k];
	    grad_te[center_l][0] += dddx;

	    dddy = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dddy += ABCD[10][k]*FourInd[k];
	    grad_te[center_l][1] += dddy;

	    dddz = 0.0;
	    for(k=0;k<quartet_size;k++)
	      dddz += ABCD[11][k]*FourInd[k];
	    grad_te[center_l][2] += dddz;

	    grad_te[center_k][0] -= ddax + ddbx + dddx;
	    grad_te[center_k][1] -= dday + ddby + dddy;
	    grad_te[center_k][2] -= ddaz + ddbz + dddz;
	}

  if (UserOptions.print_lvl >= PRINT_TEDERIV)
    print_atomvec("Two-electron contribution to the forces (a.u.)",grad_te);

  for(i=0;i<Molecule.num_atoms;i++) {
    Grad[i][0] += grad_te[i][0];
    Grad[i][1] += grad_te[i][1];
    Grad[i][2] += grad_te[i][2];
  }
  
  /*---------
    Clean-up
   ---------*/
  free(FourInd);
  free(int_stack);
  free(zero_stack);
  free(Shell_Data);
  free_block(grad_te);

  if (strcmp(UserOptions.wfn,"SCF"))
    iwl_buf_close(&TPDM,1);

  return;
}



