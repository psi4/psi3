#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <qt.h>
#include <iwl.h>
#include <math.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "defines.h"

static void make_arrays(double *evals, double **evals_docc, double **evals_virt,
                  int ndocc, int ndocc_act, int nvirt, int **ioff3, int **focact,
                  int **locact, int **fvract, int **lvract, int **vir_Q2P, int **act2full);


void mp2r12_energy(void)
{

  struct iwlbuf ERIBuff;                /* holds ERIs */
  struct iwlbuf R12Buff;                /* holds r12 integrals */
  struct iwlbuf R12TBuff;                /* holds [T1+T2,r12] integrals */

  int ndocc, ndocc_act, nvirt, nmo, ndg, nte;
  int *ioff3;
  int *clsdpi, *frdocc, *orbspi;
  int nirreps;
  int i,j,a,b,p,q,g,h;
  int ij,kl,mn,ab,pq,ip,jp,jq,iq;
  int ipjq,iqjp;
  int ij_act;
  int psym,isym,qsym,jsym,ijsym;
  int asym,bsym;
  int ifirst,jfirst,afirst,bfirst,pfirst,qfirst;
  int ilast,jlast,alast,blast,plast,qlast;
  int ntri_docc, ntri_docc_act, ntri;
  int spin;
  int print_lvl,keep_integrals;
  int *first, *last;
  int *focact,*locact,*fvract,*lvract;
  int *vir_Q2P;                         /* QTS to Pitzer index for virtual orbitals */
  int *act2full;                        /* maps active "ij" index to full "ij" index */

  double mp2_energy = 0.0;
  double mp2r12_energy = 0.0;
  double escf;
  double tolerance;
  double *eri, *r12, *r12t;
  double **K, **R, **U, **V_full[2], **V[2], **T_full[2], **B[2];
  double **Binv[2];
  double *pair_energy[2];               /* MP2 pair energies for each spin case */
  double *evals_docc, *evals_virt;
  double *tmp_ptr,**evecs,**tmp1,**tmp2;
  double spin_pfac;
  double perm_pfac;
  double VoBtrace[2], traceV, traceB;
  const double oosqrt2 = 1.0/sqrt(2.0);
  
  escf = moinfo.escf;
  nmo = moinfo.nmo;
  tolerance = params.tolerance;
  nirreps = moinfo.nirreps;
  clsdpi = moinfo.clsdpi;
  orbspi = moinfo.orbspi;
  frdocc = moinfo.frdocc;
  first = moinfo.first;
  last = moinfo.last;
  print_lvl = params.print_lvl;
  keep_integrals = params.keep_integrals;

  ndocc = 0;
  ndocc_act = 0;
  nvirt = 0;
  for(h=0; h < nirreps; h++) {
      ndocc += clsdpi[h];
      ndocc_act += clsdpi[h] - frdocc[h];
      nvirt += orbspi[h] - clsdpi[h];
    }

  ntri_docc = ndocc*(ndocc+1)/2;
  ntri_docc_act = ndocc_act*(ndocc_act+1)/2;
  ntri = nmo*(nmo+1)/2;
  ndg = ndocc*nmo;
  nte = ndg*(ndg+1)/2;
  eri = init_array(nte);
  r12 = init_array(nte);
  r12t = init_array(nte);
  pair_energy[0] = init_array(ntri_docc);
  pair_energy[1] = init_array(ntri_docc);
  K   = block_matrix(ntri_docc,ntri);
  R   = block_matrix(ntri_docc,ntri);
  U   = block_matrix(ntri_docc,ntri);
  V_full[0]   = block_matrix(ntri_docc,ntri_docc);
  V_full[1]   = block_matrix(ntri_docc,ntri_docc);
  if (ndocc != ndocc_act) {
    V[0]   = block_matrix(ntri_docc_act,ntri_docc_act);
    V[1]   = block_matrix(ntri_docc,ntri_docc);
  }
  T_full[0]   = block_matrix(ntri_docc,ntri_docc);
  T_full[1]   = block_matrix(ntri_docc,ntri_docc);
  B[0]   = block_matrix(ntri_docc_act,ntri_docc_act);
  B[1]   = block_matrix(ntri_docc_act,ntri_docc_act);
  Binv[0] = block_matrix(ntri_docc_act,ntri_docc_act);
  Binv[1] = block_matrix(ntri_docc_act,ntri_docc_act);
  tmp_ptr = init_array(ntri_docc_act);
  evecs = block_matrix(ntri_docc_act,ntri_docc_act);
  tmp1 = block_matrix(ntri_docc_act,ntri_docc_act);
  tmp2 = block_matrix(ntri_docc_act,ntri_docc_act);
  

  make_arrays(moinfo.evals,&evals_docc,&evals_virt,ndocc, ndocc_act, nmo, &ioff3,
              &focact, &locact, &fvract, &lvract, &vir_Q2P, &act2full);
  
  iwl_buf_init(&ERIBuff, PSIF_MO_TEI, tolerance, 1, 1);
  iwl_buf_rd_all_mp2r12a(&ERIBuff, eri, ioff3, ioff3, 1, ioff, (print_lvl > 3), outfile);
  iwl_buf_close(&ERIBuff, keep_integrals);

  iwl_buf_init(&R12Buff, PSIF_MO_R12, tolerance, 1, 1);
  iwl_buf_rd_all_mp2r12a(&R12Buff, r12, ioff3, ioff3, 1, ioff, (print_lvl > 3), outfile);
  iwl_buf_close(&R12Buff, keep_integrals);

  iwl_buf_init(&R12TBuff, PSIF_MO_R12T1, tolerance, 1, 1);
  iwl_buf_rd_all_mp2r12a(&R12TBuff, r12t, ioff3, ioff3, 0, ioff, (print_lvl > 3), outfile);
  iwl_buf_close(&R12TBuff, keep_integrals);

  /*------------------------------------------
    For each spin case form K, R, T matrices,
    and compute pair energies
   ------------------------------------------*/
  for(spin=0; spin <= 1; spin++) {
    spin_pfac = 1.0 - 2.0*spin;

    /*-----------------------------
      Compute K and R matrices
     -----------------------------*/
    /* NOTE: for each spin case I use only one scratch
       array for K and R since ijpq loops run over the
       same range
       */
    for(isym=0; isym < nirreps; isym++) {
      ifirst = focact[isym];
      ilast = locact[isym];
      for(i=ifirst; i <= ilast; i++) {
	for(jsym=0; jsym <= isym; jsym++) {
	  ijsym = isym^jsym;
	  jfirst = focact[jsym];
	  /*--- i >= j ---*/
	  jlast = MIN(locact[jsym],i);
	  for(j=jfirst; j <= jlast; j++) {
	    ij = ioff[i]+j;
	    for(psym=0; psym < nirreps; psym++) {
	      pfirst = first[psym];
	      plast = last[psym];
	      for(p=pfirst; p <= plast; p++) {
		ip = ioff3[i] + p;
		jp = ioff3[j] + p;
		/*--- We need p >= q -> if qsym > psym - continue ---*/
		qsym = ijsym^psym;
		if (qsym > psym)
		  continue;
		qfirst = first[qsym];
		/*--- p >= q ---*/
		qlast = MIN(last[qsym],p);
		for(q=qfirst; q <= qlast; q++) {
		  iq = ioff3[i] + q;
		  jq = ioff3[j] + q;
		  pq = ioff[p]+q;
		  ipjq = INDEX(ip,jq);
		  iqjp = INDEX(iq,jp);
		  perm_pfac = (p == q) ? oosqrt2 : 1.0;
		  perm_pfac *= (i == j) ? oosqrt2 : 1.0;
		  K[ij][pq] = perm_pfac*(eri[ipjq] + spin_pfac*eri[iqjp]);
		  R[ij][pq] = perm_pfac*(r12[ipjq] + spin_pfac*r12[iqjp]);
		  U[ij][pq] = perm_pfac*(r12t[ipjq] + spin_pfac*r12t[iqjp]);
		}
	      }
	    }
	  }
	}
      }
    }

    /*-------------------------
      Compute V and T matrices
     -------------------------*/
    for(ij=0;ij<ntri_docc;ij++) {
      V_full[spin][ij][ij] = 1.0;
      T_full[spin][ij][ij] = 1.0;
    }
    C_DGEMM('n','t',ntri_docc,ntri_docc,ntri,-1.0,
	    R[0],ntri,K[0],ntri,1.0,V_full[spin][0],ntri_docc);
    C_DGEMM('n','t',ntri_docc,ntri_docc,ntri,+1.0,          /* <- "+"-sign here because apparently I am missing it somewhere */
	    R[0],ntri,U[0],ntri,1.0,T_full[spin][0],ntri_docc);

    if (ndocc != ndocc_act) {
      for(ij=0;ij<ntri_docc_act;ij++)
	for(kl=0;kl<ntri_docc_act;kl++)
	  V[spin][ij][kl] = V_full[spin][act2full[ij]][act2full[kl]];
    }
    else
      V[spin] = V_full[spin];
    
    /*-----------------
      Compute B matrix
     -----------------*/
    for(ij=0;ij<ntri_docc_act;ij++)
      for(kl=0;kl<=ij;kl++)
	B[spin][kl][ij] = B[spin][ij][kl] = 0.5*(T_full[spin][act2full[ij]][act2full[kl]]
						 + T_full[spin][act2full[kl]][act2full[ij]]);

    free_block(T_full[spin]);
    
    /*-------------------------------
      Compute basis set completeness
     -------------------------------*/
    traceV = traceB = 0.0;
    for(i=0;i<ndocc_act;i++)
      for(j=0;j<=i-spin;j++) {
	ij = ioff[i]+j;
	traceV += V[spin][ij][ij];
	traceB += B[spin][ij][ij];
      }
    VoBtrace[spin] = traceV/traceB;

    /*----------------------
      Compute B^(-1) matrix
     ----------------------*/
    sq_rsp(ntri_docc_act,ntri_docc_act,B[spin],tmp_ptr,1,evecs,1e-14);
    for(i=0;i<ntri_docc_act;i++)
      tmp1[i][i] = 1.0/tmp_ptr[i];
    mmult(tmp1,0,evecs,1,tmp2,0,ntri_docc_act,ntri_docc_act,ntri_docc_act,0);  
    mmult(evecs,0,tmp2,0,Binv[spin],0,ntri_docc_act,ntri_docc_act,ntri_docc_act,0);
    
    /*---------------------------------------
      Compute conventional MP2 pair energies
     ---------------------------------------*/
    spin_pfac = 2.0*spin + 1.0;
    for(isym=0; isym < nirreps; isym++) {
      ifirst = focact[isym];
      ilast = locact[isym];
      for(i=ifirst; i <= ilast; i++) {
	for(jsym=0; jsym <= isym; jsym++) {
	  jfirst = focact[jsym];
	  /*--- i >= j ---*/
	  jlast = MIN(locact[jsym],i);
	  ijsym = isym^jsym;
	  for(j=jfirst; j <= jlast; j++) {
	    ij = INDEX(i,j);
	    for(asym=0; asym < nirreps; asym++) {
	      afirst = fvract[asym];
	      alast = lvract[asym];
	      for(a=afirst; a <= alast; a++) {
		/*--- We need a >= b -> if bsym > asym - continue ---*/
		bsym = ijsym^asym;
		if (bsym > asym)
		  continue;
		bfirst = fvract[bsym];
		blast = MIN(lvract[bsym],a);
		for(b=bfirst; b <= blast; b++) {
		  ab = INDEX(vir_Q2P[a],vir_Q2P[b]);
		  pair_energy[spin][ij] += spin_pfac*K[ij][ab]*K[ij][ab]/
			    (evals_docc[i] + evals_docc[j] -
			     evals_virt[a] - evals_virt[b]);
		}
	      }
	    }
	  }
	}
      }
    }

    if (print_lvl > 1) {
      fprintf(outfile," %s V-matrix (the difference between the unit matrix and its\n r12/r12 representation in this basis):\n",
	      (spin == 0) ? "Singlet" : "Triplet");
      print_mat(V[spin],ntri_docc_act,ntri_docc_act,outfile);
      fprintf(outfile,"\n");

      fprintf(outfile," %s B-matrix (the difference between the unit matrix and its\n ???? representation in this basis)::\n",
	      (spin == 0) ? "Singlet" : "Triplet");
      print_mat(B[spin],ntri_docc_act,ntri_docc_act,outfile);
      fprintf(outfile,"\n");
    }

    for(ij=0;ij<ntri_docc_act;ij++)
      mp2_energy += pair_energy[spin][act2full[ij]];


    /*-------------------------------------------------
      Compute MP2-R12/A contributions to pair energies
     -------------------------------------------------*/
    spin_pfac = 2.0*spin + 1.0;
    for(i=0;i<ndocc_act;i++)
      for(j=0;j<=i-spin;j++) {
	ij = ioff[i]+j;
	for(kl=0;kl<ntri_docc_act;kl++)
	  for(mn=0;mn<ntri_docc_act;mn++)
	    if (!params.c_limit)
	      pair_energy[spin][act2full[ij]] -= spin_pfac*V[spin][kl][ij]*V[spin][mn][ij]*Binv[spin][kl][mn];
	    else
	      pair_energy[spin][act2full[ij]] -= spin_pfac*V[spin][ij][ij]*(spin == 0 ? -0.5 : -0.25);
      }
    for(ij=0;ij<ntri_docc_act;ij++)
      mp2r12_energy += pair_energy[spin][act2full[ij]];

    fprintf(outfile,"\t%s pair energies:\n",(spin == 0) ? "Singlet" : "Triplet");
    fprintf(outfile,"\t    i       j         e(ij)\n");
    fprintf(outfile,"\t  -----   -----   ------------\n");
    ij_act = 0;
    for(g=0; g < nirreps; g++)
    for(i=focact[g]; i <= locact[g] ; i++) {
      for(h=0; h <= g; h++) {
	jlast  = (locact[h] > i) ? i : locact[h];
	for(j=focact[h]; j <= jlast ; j++,ij_act++) {
	  if (spin == 0 || i != j)
	    fprintf(outfile,"\t  %3d     %3d     %12.9lf\n",i+1,j+1,pair_energy[spin][act2full[ij_act]]);
	}
      }
    }
    fprintf(outfile,"\n");

    fprintf(outfile,"\t%s -Tr(V)/Tr(B) = %lf\n\n",(spin == 0) ? "Singlet" : "Triplet",
	    (-1.0)*VoBtrace[spin]);
    
    
  }

  fprintf(outfile, "\n\tMBPT(2) Energy        = %20.10lf\n", mp2_energy);
  fprintf(outfile,   "\tMBPT(2)-R12/A Energy  = %20.10lf\n", mp2r12_energy);
  fprintf(outfile,   "\tTotal Energy          = %20.10lf\n", escf+mp2r12_energy);
  fflush(outfile);

  free(eri);
  free(r12);
  free_block(K);
  free_block(R);
  if (ndocc_act != ndocc) {
    free_block(V[0]);
    free_block(V[1]);
  }
  free(pair_energy[0]);
  free(pair_energy[1]);
  free(tmp_ptr);
  free_block(evecs);
  free_block(tmp1);
  free_block(tmp2);

}

void make_arrays(double *evals, double **evals_docc, double **evals_virt,
                 int ndocc, int ndocc_act, int nmo, int **ioff3, int **focact,
                 int **locact, int **fvract, int **lvract, int **vir_Q2P, int **act2full)
{

  int g,h,i,j,iact,jact;
  int jlast;
  int num_docc;
  int count, virt_count, offset;
  int *clsdpi;
  int *orbspi;
  int *virtpi;
  int *frdocc;
  int *fruocc;
  int nirreps;
  int first_offset, last_offset;

  clsdpi = moinfo.clsdpi;
  orbspi = moinfo.orbspi;
  nirreps = moinfo.nirreps;
  virtpi = moinfo.virtpi;
  frdocc = moinfo.frdocc;
  fruocc = moinfo.fruocc;

  /* Construct first and last indexing arrays for occupied and virtual
     orbitals in QTP ordering */

  *focact = init_int_array(nirreps);
  *locact = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
      (*focact)[h] = -1;
      (*locact)[h] = -2;
    }
  first_offset=frdocc[0];
  last_offset=clsdpi[0] - 1;
  (*focact)[0] = first_offset;
  (*locact)[0] = last_offset;
  for(h=1; h < nirreps; h++) {
      first_offset += clsdpi[h-1] - frdocc[h-1] + frdocc[h];
      last_offset += clsdpi[h];
      if(clsdpi[h]-frdocc[h]) {
          (*focact)[h] = first_offset;
          (*locact)[h] = last_offset;
        }
    }

  *fvract = init_int_array(nirreps);
  *lvract = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
      (*fvract)[h] = -1;
      (*lvract)[h] = -2;
    }
  first_offset=0;
  last_offset=virtpi[0] - fruocc[0] - 1;
  (*fvract)[0] = first_offset;
  (*lvract)[0] = last_offset;
  for(h=1; h < nirreps; h++) {
      first_offset += virtpi[h-1];
      last_offset += virtpi[h] - fruocc[h] + fruocc[h-1];
      if(virtpi[h] - fruocc[h]) {
          (*fvract)[h] = first_offset;
          (*lvract)[h] = last_offset;
        }
    }

  (*vir_Q2P) = init_int_array(nmo-ndocc);
  virt_count = 0;
  offset = 0;
  for(h=0; h < nirreps; h++) {
    if (h)
      offset += orbspi[h-1];
    for(i=clsdpi[h]; i < orbspi[h]; i++) {
      (*vir_Q2P)[virt_count] = i + offset;
      virt_count++;
    }
  }

  /* Re-order eigenvalues from Pitzer to QTP Relative ordering */
  (*evals_docc) = init_array(ndocc);
  count=0;
  offset=0;
  for(h=0; h < nirreps; h++) {
      if(h) offset += orbspi[h-1];
      num_docc = clsdpi[h];
      for(j=offset; j < (offset + num_docc); j++) {
          (*evals_docc)[count] = evals[j];
          count++;
        }
    }

  (*evals_virt) = init_array(nmo-ndocc);
  count=0;
  offset=0;
  for(h=0; h < nirreps; h++) {
      if(h) offset += orbspi[h-1];
      for(j=offset+clsdpi[h]; j < (offset+orbspi[h]); j++) {
          (*evals_virt)[count] = evals[j];
          count++;
        }
    }

  /* Generate ioff3 array.  This array gives the row offset for an
     ndocc x nmo matrix */
  (*ioff3) = init_int_array(MAXIOFF3);
  for(i=0; i < MAXIOFF3; i++) {
      (*ioff3)[i] = i*nmo;
    }

  /* Generate active "ij" to full "ij" mapping */
  (*act2full) = init_int_array(ndocc_act*(ndocc_act+1)/2);
  iact = 0;
  for(g=0; g < nirreps; g++)
    for(i=(*focact)[g]; i <= (*locact)[g] ; i++,iact++) {
      jact = 0;
      for(h=0; h <= g; h++) {
	jlast  = ((*locact)[h] > i) ? i : (*locact)[h];
	for(j=(*focact)[h]; j <= jlast ; j++,jact++) {
	  (*act2full)[ioff[iact]+jact] = ioff[i] + j;
	}
      }
    }
  
  return;
}


