#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <ip_libv1.h>
#include <qt.h>
#include <iwl.h>
#include <file30.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"


#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void trans_one_forwards(void);
void trans_one_backwards(void);
void tran_one(int nirreps, double ***C, int src_orbs, int *src_first, int *src_last,
              int *src_orbspi, int dst_orbs, int *dst_first, int *dst_last,
              int *dst_orbspi, double *src_ints, double *dst_ints, int *order,
              char *label, int backtran, int nfzc, int printflg,FILE *outfile);
double *** construct_evects(char *spin, int nirreps, int *active, int *sopi, int *orbspi,
                            int *first_so, int *last_so, int *first, int *last,
			    int *fstact, int *lstact, int printflag);

/*
** TRANSFORM_ONE(): This function does the SO to MO transformation of
**    the one-electron integrals.
**
** This function was getting too complicated to read, so I split
**  it up into two functions, one for forwards transforms, and one
**  for backwards transforms.
**
** C. David Sherrill
** 29 April 1998
**
*/
void transform_one()
{

  if (!strcmp(params.wfn,"MP2")) return;
  if (params.backtr) trans_one_backwards();
  else trans_one_forwards();

}



/*
** TRANS_ONE_FORWARDS()
**
** This does a forwards transform of one-electron integrals from the SO
**  to the MO basis.  All the real work is done in a function
**  tran_one(); the rest is just setting up to call that function. 
**  I also simplified matters by changing libiwl to allow filtering
**  out of frozen orbitals.  Hence, for forwards transforms, we now
**  only write (up to) two files: the bare one-electron integrals,
**  and the frozen core operator.  Each file runs over all indices 
**  (frozen or not).  This makes this code much simpler than it was.
**
** C. David Sherrill
** 29 April 1998
*/
void trans_one_forwards(void)
{
  int itap, nirreps, nmo, nso, nfzc, nfzv;
  int *src_first, *dst_first, *src_last, *dst_last, *src_orbspi, *dst_orbspi;
  int src_orbs, dst_orbs, src_ntri, dst_ntri;
  int print_integrals;
  double *oe_ints, *new_oe_ints;

  if (params.print_lvl) {
    fprintf(outfile, "\n\tTransforming one-electron integrals...\n");
    fflush(outfile);
  }

  nirreps = moinfo.nirreps;
  nmo = moinfo.nmo;
  nso = moinfo.nso;
  print_integrals = params.print_oe_ints;
  nfzc = 0;  nfzv = 0;

  src_first = moinfo.first_so;
  dst_first = moinfo.first;
  src_last = moinfo.last_so;
  dst_last = moinfo.last;
  src_orbspi = moinfo.sopi;
  dst_orbspi = moinfo.orbspi;

  src_orbs = nso;
  dst_orbs = nmo - nfzv - nfzc;
  dst_ntri = (dst_orbs * (dst_orbs + 1)) / 2;
  src_ntri = (src_orbs * (src_orbs + 1)) / 2;

  new_oe_ints = init_array(dst_ntri);

  /* We do the two-electron ints first, for which our target ints may
   * not have allowed frozen indices.  However, here we now want to run
   * over all MOs.  So, we need to re-build the C matrix to include 
   * frozen orbital columns.  This annoying stuff could be avoided by
   * a generalized mmult routine
   */ 
  if (!params.do_all_tei) {
    destruct_evects(nirreps, moinfo.evects);
    destruct_evects(nirreps, moinfo.evects_alpha);
    destruct_evects(nirreps, moinfo.evects_beta);

    moinfo.evects = construct_evects("alpha", moinfo.nirreps, moinfo.orbspi,
				     moinfo.sopi, moinfo.orbspi,
				     moinfo.first_so, moinfo.last_so,
				     moinfo.first, moinfo.last,
				     moinfo.first, moinfo.last, params.print_mos);
    moinfo.evects_alpha = construct_evects("alpha", moinfo.nirreps, moinfo.orbspi,
					   moinfo.sopi, moinfo.orbspi,
					   moinfo.first_so, moinfo.last_so,
					   moinfo.first, moinfo.last,
					   moinfo.first, moinfo.last, params.print_mos);
    moinfo.evects_beta = construct_evects("beta", moinfo.nirreps, moinfo.orbspi,
					  moinfo.sopi, moinfo.orbspi,
					  moinfo.first_so, moinfo.last_so,
					  moinfo.first, moinfo.last,
					  moinfo.first, moinfo.last, params.print_mos);
  }

  if (params.do_h_bare) {
    
    oe_ints = moinfo.oe_ints;

    if(!strcmp(params.ref,"UHF")) {
    
      itap = params.h_bare_a_file;

      tran_one(nirreps, moinfo.evects_alpha, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	       dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order_alpha,
	       "\n\n\tOne-Electron Alpha Integrals (MO Basis):\n", 
	       params.backtr, nfzc, print_integrals, outfile);

      iwl_wrtone(itap, dst_ntri, new_oe_ints, moinfo.efzc);

      if (params.print_lvl) {
	fprintf(outfile, "\tOne-electron A integrals written to file %d.\n", itap);
	fflush(outfile);
      }
      itap = params.h_bare_b_file;

      /* Clear the new_oe_ints */
      zero_arr(new_oe_ints,dst_ntri);

      tran_one(nirreps, moinfo.evects_beta, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	       dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order_beta,
	       "\n\n\tOne-Electron Beta Integrals (MO Basis):\n", 
	       params.backtr, nfzc, print_integrals, outfile);

      iwl_wrtone(itap, dst_ntri, new_oe_ints, moinfo.efzc);
      if (params.print_lvl) {
	fprintf(outfile, "\tOne-electron B integrals written to file %d.\n", itap);
	fflush(outfile);
      }
    }
    else {

      itap = params.h_bare_file;

      tran_one(nirreps, moinfo.evects, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	       dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order,
	       "\n\n\tOne-Electron Integrals (MO Basis):\n", 
	       params.backtr, nfzc, print_integrals, outfile);


      iwl_wrtone(itap, dst_ntri, new_oe_ints, moinfo.efzc);
      if (params.print_lvl) {
	fprintf(outfile, "\tOne-electron integrals written to file %d.\n", itap);
	fflush(outfile);
      }
    }


  }


  if (params.do_h_fzc) {

    if(!strcmp(params.ref,"UHF")) {

      itap = params.h_fzc_a_file;
      oe_ints = moinfo.fzc_operator_alpha;

      tran_one(nirreps, moinfo.evects_alpha, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	       dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order_alpha,
	       "\n\n\tAlpha Frozen-Core Operator (MO Basis):\n", 
	       params.backtr, nfzc, print_integrals, outfile);

      iwl_wrtone(itap, dst_ntri, new_oe_ints, moinfo.efzc);

      if (params.print_lvl) {
	fprintf(outfile, "\tAlpha frozen-core operator written to file %d.\n", itap);
	fflush(outfile);
      }

      /* Clear the new_oe_ints */
      zero_arr(new_oe_ints,dst_ntri);

      itap = params.h_fzc_b_file;
      oe_ints = moinfo.fzc_operator_beta;

      tran_one(nirreps, moinfo.evects_beta, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	       dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order_beta,
	       "\n\n\tBeta Frozen-Core Operator (MO Basis):\n", 
	       params.backtr, nfzc, print_integrals, outfile);

      iwl_wrtone(itap, dst_ntri, new_oe_ints, moinfo.efzc);

      if (params.print_lvl) {
	fprintf(outfile, "\tBeta frozen-core operator written to file %d.\n", itap);
	fflush(outfile);
      }
    }
    else {

      itap = params.h_fzc_file;
      oe_ints = moinfo.fzc_operator;

      tran_one(nirreps, moinfo.evects, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	       dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order,
	       "\n\n\tFrozen-Core Operator (MO Basis):\n", 
	       params.backtr, nfzc, print_integrals, outfile);

      iwl_wrtone(itap, dst_ntri, new_oe_ints, moinfo.efzc);

      if (params.print_lvl) {
	fprintf(outfile, "\tFrozen-core operator written to file %d.\n", itap);
	fflush(outfile);
      }

    }

  }

  free(new_oe_ints);
  
}
 

/*
** TRANS_ONE_BACKWARDS()
**
** This does a backwards transform of the one-particle density matrix
**  and the MO Lagrangian from the MO to the AO basis.  All the real 
**  work is done in a function tran_one()---the same one called in the
**  forwards transformations; the code in this function just sets up 
**  to call tran_one().
**
** C. David Sherrill
** 29 April 1998
*/
void trans_one_backwards(void)
{
  int itap, nirreps, nmo, nfzc, nfzv;
  int *src_first, *dst_first, *src_last, *dst_last, *src_orbspi, *dst_orbspi;
  int src_orbs, dst_orbs, src_ntri, dst_ntri;
  int print_integrals;
  double *oe_ints, *new_oe_ints;
  double **A,**B,*opdm,**tmat,**tmat2,**so2ao,tval;
  int p,q,pq;
  int I,J,P,Q,PQ;
  PSI_FPTR pdm_fp_in=0, pdm_fp_out=0;

  if (params.print_lvl) {
    fprintf(outfile, "\n\tTransforming one-electron integrals...\n");
    fflush(outfile);
  }

  itap = params.opdm_out_file;
  nirreps = moinfo.backtr_nirreps;
  nmo = moinfo.nmo;
  print_integrals = params.print_oe_ints;
  nfzc = 0;
  nfzv = moinfo.nfzv;   
  
  src_first = moinfo.backtr_mo_first;
  dst_first = moinfo.backtr_ao_first;
  src_last = moinfo.backtr_mo_lstact;
  dst_last = moinfo.backtr_ao_last;
  src_orbspi = moinfo.backtr_mo_active;
  dst_orbspi = moinfo.backtr_ao_orbspi;
  
  src_orbs = nmo - nfzv;
  dst_orbs = moinfo.nao;
  dst_ntri = (dst_orbs * (dst_orbs + 1)) / 2;
  src_ntri = (src_orbs * (src_orbs + 1)) / 2;

  /* read in the MO one-particle density matrix */
  opdm = init_array(src_ntri);

  tmat = block_matrix(src_orbs, src_orbs);

  rfile(params.opdm_in_file);
  if (flen(params.opdm_in_file)==0) {
    fprintf(outfile,"\tCan't find the one-particle density matrix file %d\n",
	    params.opdm_in_file);
    rclose(params.opdm_in_file, 4);
    return;
  }
  if (strcmp(params.wfn, "OOCCD")!=0)
    wreadw(params.opdm_in_file, (char *) tmat[0], src_orbs * src_orbs * 
           sizeof(double), pdm_fp_in, &pdm_fp_in);
  else {
    tmat2 = block_matrix(nmo, nmo);
    wreadw(params.opdm_in_file, (char *) tmat2[0], nmo * nmo *
           sizeof(double), pdm_fp_in, &pdm_fp_in);
    for (p=0; p<src_orbs; p++) 
      for (q=0; q<src_orbs; q++) 
        tmat[p][q] = tmat2[p][q];
    free_block(tmat2);
  } 

  if (params.print_lvl > 3) {
    fprintf(outfile, "One-particle density matrix\n");
    print_mat(tmat, src_orbs, src_orbs, outfile);
  }
  rclose(params.opdm_in_file, 3);
  
  /* reorder the opdm into the backtransform order */
  for (p=0; p<src_orbs; p++) {
    for (q=0; q<=p; q++) {
      P = moinfo.corr2pitz_nofzv[p];
      Q = moinfo.corr2pitz_nofzv[q];
      PQ = INDEX(P,Q);
      opdm[PQ] = 0.5 * (tmat[p][q] + tmat[q][p]);
    }
  } 
  
  oe_ints = opdm;
  new_oe_ints = init_array(dst_ntri);
  
  tran_one(nirreps, moinfo.evects, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
	   dst_first, dst_last, dst_orbspi, oe_ints, new_oe_ints, moinfo.order,
           "\n\n\tOne-Particle Density Matrix (AO Basis):\n", 
           params.backtr, nfzc, print_integrals, outfile);

  rfile(itap);
  wwritw(itap, (char *) new_oe_ints, dst_ntri*sizeof(double), pdm_fp_out,
	 &pdm_fp_out);
  free(opdm);
  free_block(tmat);


  /* If this is a back-transform, we need to do the Lagrangian also */
  /* Assume the lagrangian has indices over all orbitals, incl fzv  */
  src_last = moinfo.backtr_mo_last;
  src_orbspi = moinfo.backtr_mo_orbspi;
  src_orbs = nmo; 
  src_ntri = (src_orbs * (src_orbs + 1)) / 2;
  
  /* we need to re-construct the C matrix because we now need the   */
  /* columns corresponding to frozen virtual orbitals               */
  destruct_evects(moinfo.backtr_nirreps, moinfo.evects);
  moinfo.evects = (double ***) malloc (1 * sizeof(double **));
  moinfo.evects[0] = block_matrix(moinfo.nao, moinfo.nmo);
  file30_init();
  so2ao = file30_rd_usotao_new();
  file30_close();
  if (params.print_mos) print_mat(so2ao,moinfo.nso,moinfo.nao,outfile);
  mmult(so2ao,1,moinfo.scf_vector,0,moinfo.evects[0],0,moinfo.nao,moinfo.nso,
	moinfo.nmo,0);
  if (params.print_mos) {
    fprintf(outfile, "C matrix (including AO to SO)\n");
    print_mat(moinfo.evects[0], moinfo.nao, moinfo.nmo, outfile);
  }
/* TDC --- converted libfile30 to use block_matrix() */
  free_block(so2ao);
  tmat = block_matrix(src_orbs, src_orbs);


  /* read in the MO Lagrangian */
  opdm = init_array(src_ntri);
  rfile(params.lag_in_file);
  if (flen(params.lag_in_file)==0) {
    fprintf(outfile,"\tCan't find the Lagrangian matrix file %d\n",
	    params.lag_in_file);
    rclose(params.lag_in_file, 4);
    return;
  }
  pdm_fp_in = 0;
  wreadw(params.lag_in_file,(char *) tmat[0],
	 src_orbs*src_orbs*sizeof(double), pdm_fp_in, &pdm_fp_in);
  if (params.print_lvl > 3) {
    fprintf(outfile, "Lagrangian (MO Basis):\n");
    print_mat(tmat, src_orbs, src_orbs, outfile);
  }
  rclose(params.lag_in_file, 3);
  
  /* reorder the Lagrangian to the back-transform order 
   * and enforce permutational symmetry.  
   */
  for (p=0; p<src_orbs; p++) {
    for (q=0; q<=p; q++) {
      P = moinfo.corr2pitz[p];
      Q = moinfo.corr2pitz[q];
      tval = 0.5 * (tmat[p][q] + tmat[q][p]);
      if (params.lagran_double) tval *= 2.0;
      if (params.lagran_halve) tval *= 0.5;
      PQ = INDEX(P,Q);
      opdm[PQ] += tval;
    }
  } 
  free_block(tmat);
  
  if (params.print_lvl > 3) {
    fprintf(outfile, "Reordered, Symmetrized Lagrangian in MO basis\n");
    print_array(opdm, src_orbs, outfile);
  }
  
  /* now do the backtransformation */
  tran_one(nirreps, moinfo.evects, src_orbs, src_first, src_last, src_orbspi, dst_orbs, 
           dst_first, dst_last, dst_orbspi, opdm, new_oe_ints, moinfo.order,
           "\n\n\tLagrangian (AO Basis):\n", 
	   params.backtr, nfzc, print_integrals, outfile);
  
  /* append the AO Lagrangian to the AO one-particle density matrix */
  wwritw(itap, (char *) new_oe_ints, dst_ntri*sizeof(double), pdm_fp_out,
	 &pdm_fp_out);
  rclose(itap,3);
  free(opdm);


  if (params.print_lvl) {
    fprintf(outfile, "\tOne-pdm and lagrangian written to file%d.\n", itap);
  }

  free(new_oe_ints);
  
}
 

 

void tran_one(int nirreps, double ***C, int src_orbs, int *src_first, int *src_last,
              int *src_orbspi, int dst_orbs, int *dst_first, int *dst_last,
              int *dst_orbspi, double *src_ints, double *dst_ints, int *order,
              char *label, int backtran, int nfzc, int printflg, FILE *outfile)
{

  int psym, p, q, P, Q, pfirst, plast, pq;
  int i, j, ifirst, ilast, I, J, i2, j2, ij;
  double **A, **B;
  int A_cols, B_cols, *C_cols;
  
  A_cols = MAX0(src_orbs,dst_orbs);
  A = block_matrix(A_cols, A_cols);
  B = block_matrix(src_orbs,dst_orbs);
  B_cols = dst_orbs;
  C_cols = backtran ? src_orbspi : dst_orbspi;

  if (printflg) {
    fprintf(outfile, "%s\n", label);
  }

  for (psym=0; psym < nirreps; psym++) {
    pfirst = src_first[psym];
    plast = src_last[psym];
    for (p=pfirst,P=0; p <= plast; p++,P++) {
      for (q=pfirst,Q=0; q <= plast; q++,Q++) {
        pq = INDEX(p,q);
        A[P][Q] = src_ints[pq];
      }
    }
    
    if (src_orbspi[psym] && dst_orbspi[psym]) {

#ifdef USE_BLAS
      C_DGEMM('n', backtran ? 't' : 'n', src_orbspi[psym], dst_orbspi[psym],
              src_orbspi[psym], 1.0, A[0], A_cols, C[psym][0], C_cols[psym],
              0.0, B[0], B_cols);

      C_DGEMM(backtran ? 'n' : 't', 'n', dst_orbspi[psym], dst_orbspi[psym],
              src_orbspi[psym], 1.0, C[psym][0], C_cols[psym], B[0], B_cols,
              0.0, A[0], A_cols);
#else
      mmult(A,0,C[psym],(backtran ? 1 : 0),B,0,
            src_orbspi[psym],src_orbspi[psym],dst_orbspi[psym],0);
      mmult(C[psym],(backtran ? 0 : 1),B,0,A,0,
            dst_orbspi[psym],src_orbspi[psym],dst_orbspi[psym],0);
#endif

    }    

    if (printflg) {
      if (dst_orbspi[psym]) {
        fprintf(outfile, " Irrep %s\n", moinfo.labels[psym]);
        print_mat(A,dst_orbspi[psym],dst_orbspi[psym],outfile);
      }
    }

    ifirst = dst_first[psym];
    ilast = dst_last[psym];
    for (i=ifirst,I=0; i <= ilast; i++,I++) {
      for (j=ifirst,J=0; (j <= ilast) && (j <= i); j++, J++) {

        if (!backtran) {
          i2 = order[i]-nfzc;
          j2 = order[j]-nfzc;
          ij = INDEX(i2,j2);
        }
        else
          ij = INDEX(i,j);

        dst_ints[ij] = A[I][J];
      }
    }
  } /* end loop over irreps */

  free_block(A);
  free_block(B);

}


