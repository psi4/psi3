#include <stdio.h>
#include <math.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"
#include <physconst.h>

void sigmaSS(int index, int C_irr);
void precondition_SS(dpdfile2 *RIA, dpdfile2 *Ria, double eval);
void schmidt_add_SS(dpdfile2 *RIA, dpdfile2 *Ria, int *numCs);
void precondition_SS_RHF(dpdfile2 *RIA, double eval);
void schmidt_add_SS_RHF(dpdfile2 *RIA, int *numCs);
void restart_SS(double **alpha, int L, int num, int C_irr);
void dgeev_eom(int L, double **G, double *evals, double **alpha);
double norm_C1(dpdfile2 *C1A, dpdfile2 *C1B);
double scm_C1(dpdfile2 *C1A, dpdfile2 *C1B, double a);
void init_S2(int index, int C_irr);
void init_C2(int index, int C_irr);

void diagSS(int C_irr) {
  dpdfile2 Fmi, FMI, Fae, FAE, Fme, FME;
  dpdfile2 CME, Cme, C, SIA, Sia, RIA, Ria;
  dpdbuf4 CMNEF, Cmnef, CMnEf, W;
  char lbl[32], lbl2[32];
  int lwork, info, get_right_ev = 1, get_left_ev = 0;
  int L,h,i,j,k,a,C_index,errcod,keep_going=1,numCs,iter=0;
  double norm, tval, *lambda, *lambda_old;
  double **G, *work, *evals_complex, **alpha, **evectors_left;
  int nirreps, *occpi, *openpi, *virtpi, range;
  int *converged, num_converged, num_roots;
  int begin_occ, begin_virt, end_virt, dim_SS = 0;
  int pf, cnt, irr_occ, irr_virt;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;
  range = eom_params.excitation_range;
  pf = eom_params.print_singles;

  if (pf) fprintf(outfile,"\n\n");

  /* Setup initial C guess vector(s) */
  C_index=0;

  /* This is just a bunch of tedious code to setup reasonable HOMO-LUMO
     guess vectors */

  range = eom_params.rpi[C_irr];
  for (cnt=0; cnt<nirreps; ++cnt) {
    irr_occ = dpd_dp[C_irr][cnt][0];
    irr_virt = dpd_dp[C_irr][cnt][1];

    if(params.ref == 0) { /** RHF **/

      begin_occ = occpi[irr_occ]-range;
      if (begin_occ<0) begin_occ = 0;
      end_virt = range;
      if (end_virt > virtpi[irr_virt])
	end_virt = virtpi[irr_virt];
      for (i=begin_occ; i < occpi[irr_occ]; ++i) {
	for (a=0; a < end_virt; ++a, ++C_index) {
	  sprintf(lbl, "%s %d", "CME", C_index);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&CME);
	  CME.matrix[cnt][i][a] = 1.0/sqrt(2.0);
	  dpd_file2_mat_wrt(&CME);
	  dpd_file2_mat_close(&CME);
	  dpd_file2_close(&CME);

	  sprintf(lbl, "%s %d", "Cme", C_index);
	  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&Cme);
	  Cme.matrix[cnt][i][a] = 1.0/sqrt(2.0);
	  dpd_file2_mat_wrt(&Cme);
	  dpd_file2_mat_close(&Cme);
	  dpd_file2_close(&Cme);
	}
      }

    }
    else { /** ROHF/UHF **/

      begin_occ = occpi[irr_occ]-range;
      if (begin_occ<0) begin_occ = 0;
      end_virt = range;
      if (end_virt > (virtpi[irr_virt]-openpi[irr_virt]) )
	end_virt = virtpi[irr_virt]-openpi[irr_virt];
      for (i=begin_occ; i < occpi[irr_occ]; ++i) {
	for (a=0; a < end_virt; ++a, ++C_index) {
	  sprintf(lbl, "%s %d", "CME", C_index);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&CME);
	  CME.matrix[cnt][i][a] = 1.0;
	  dpd_file2_mat_wrt(&CME);
	  dpd_file2_close(&CME);
	  sprintf(lbl, "%s %d", "Cme", C_index);
	  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&Cme);
	  dpd_file2_mat_wrt(&Cme);
	  dpd_file2_close(&Cme);
	}
      }
      begin_occ = occpi[irr_occ]-openpi[irr_occ]-range;
      if (begin_occ<0) begin_occ = 0;
      begin_virt = virtpi[irr_virt] - openpi[irr_virt];
      end_virt = virtpi[irr_virt];
      if (openpi[irr_virt] > range) end_virt = virtpi[irr_virt] - range;
      for (i=begin_occ; i < occpi[irr_occ] - openpi[irr_occ]; ++i) {
	for (a=begin_virt; a < end_virt; ++a, ++C_index) {
	  sprintf(lbl, "%s %d", "Cme", C_index);
	  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&Cme);
	  Cme.matrix[cnt][i][a] = 1.0;
	  dpd_file2_mat_wrt(&Cme);
	  dpd_file2_close(&Cme);

	  sprintf(lbl, "%s %d", "CME", C_index);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&CME);
	  dpd_file2_mat_wrt(&CME);
	  dpd_file2_close(&CME);
	}
      }
      begin_occ = occpi[irr_occ]-openpi[irr_occ]-range;
      if (begin_occ<0) begin_occ = 0;
      end_virt = range - openpi[irr_virt];
      if (end_virt<0) end_virt = 0;
      for (i=begin_occ; i < occpi[irr_occ] - openpi[irr_occ]; ++i) {
	for (a=0; a < end_virt; ++a, ++C_index) {
	  sprintf(lbl, "%s %d", "Cme", C_index);
	  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&Cme);
	  Cme.matrix[cnt][i][a] = 1.0;
	  dpd_file2_mat_wrt(&Cme);
	  dpd_file2_close(&Cme);

	  sprintf(lbl, "%s %d", "CME", C_index);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_mat_init(&CME);
	  dpd_file2_mat_wrt(&CME);
	  dpd_file2_close(&CME);
	}
      }
    }
  }

  /* Setup residual vector file */
  dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
  dpd_file2_mat_init(&RIA);
  dpd_file2_mat_wrt(&RIA);
  dpd_file2_close(&RIA);

  if(params.ref) {
    dpd_file2_init(&Ria, EOM_R, C_irr, 0, 1, "Ria");
    dpd_file2_mat_init(&Ria);
    dpd_file2_mat_wrt(&Ria);
    dpd_file2_close(&Ria);
  }

  L = num_roots = C_index;
  /* arrays must be dimensioned with at least the final number of roots - even though
     num_roots may be limited until the first collapse by the number of good
     initial guesses obtained above. */

  if (eom_params.rpi[C_irr] > C_index) i = eom_params.rpi[C_irr];
  else i = C_index;

  converged = init_int_array(i);
  lambda_old = init_array(i);

  /* make zero C2's and S2's (needed by sigmaSS's checksum) */
  for (i=0; i < (eom_params.vectors_per_root_SS+1) * eom_params.rpi[C_irr]; ++i) {
    init_C2(i, C_irr);
    init_S2(i, C_irr);
  }

  while ((keep_going == 1) && (iter < eom_params.max_iter_SS)) {
    if (pf) fprintf(outfile,"Iter=%-4d L=%-4d", iter+1, L);

    keep_going = 0;
    numCs = L;
    num_converged = 0;

    for (i=0;i<L;++i) {
      if(params.ref == 0) {
	sprintf(lbl, "%s %d", "SIA", i);
	dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	dpd_file2_scm(&SIA, 0.0);
	dpd_file2_close(&SIA);
      }
      else {
	sprintf(lbl, "%s %d", "SIA", i);
	dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	sprintf(lbl, "%s %d", "Sia", i);
	dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
	scm_C1(&SIA, &Sia, 0.0);
	dpd_file2_close(&SIA);
	dpd_file2_close(&Sia);
      }
    }

    for (i=0;i<L;++i)
      sigmaSS(i,C_irr);

    G = block_matrix(L,L);

    for (i=0;i<L;++i) {
      if(params.ref == 0) {
	sprintf(lbl, "%s %d", "CME", i);
	dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	dpd_file2_copy(&CME, CC_TMP0, "CME");
	dpd_file2_close(&CME);

	dpd_file2_init(&CME, CC_TMP0, C_irr, 0, 1, "CME");
	dpd_file2_scm(&CME, 2.0);
      }
      else {
	sprintf(lbl, "%s %d", "CME", i);
	dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	sprintf(lbl, "%s %d", "Cme", i);
	dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
      }
      for (j=0;j<L;++j) {
	if(params.ref == 0) {
	  sprintf(lbl, "%s %d", "SIA", j);
	  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	  tval = dpd_file2_dot(&CME, &SIA);
	  dpd_file2_close(&SIA);
	}
	else {
	  sprintf(lbl, "%s %d", "SIA", j);
	  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	  tval = dpd_file2_dot(&CME, &SIA);
	  dpd_file2_close(&SIA);
	  sprintf(lbl, "%s %d", "Sia", j);
	  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
	  tval += dpd_file2_dot(&Cme, &Sia);
	  dpd_file2_close(&Sia);
	}
	G[i][j] = tval;
      }
      dpd_file2_close(&CME);
      if(params.ref) dpd_file2_close(&Cme);
    }

    /*    print_mat(G, L, L, outfile); */

    lambda = init_array(L);        /* holds real part of eigenvalues of G */
    alpha = block_matrix(L,L);     /* will hold eigenvectors of G */

    dgeev_eom(L, G, lambda, alpha);

    eigsort(lambda, alpha, L);

    /*    eivout(alpha, lambda, L, L, outfile);  */

    free_block(G);

    if (pf) fprintf(outfile,
		    "  Root    EOM Energy     Delta E   Res. Norm    Conv?\n");
    fflush(outfile);

    dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
    if(params.ref) dpd_file2_init(&Ria, EOM_R, C_irr, 0, 1, "Ria");

    for (k=0;k<num_roots;++k) {
      if(params.ref == 0) dpd_file2_scm(&RIA, 0.0);
      else {
	dpd_file2_scm(&RIA, 0.0);
	dpd_file2_scm(&Ria, 0.0);
      }
      converged[k] = 0;
      for (i=0;i<L;++i) { 
	if(params.ref == 0) {
	  sprintf(lbl, "%s %d", "SIA", i);
	  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	  sprintf(lbl, "%s %d", "CME", i);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_axpbycz(&CME, &SIA, &RIA, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
	  dpd_file2_close(&CME);
	  dpd_file2_close(&SIA);
	}
	else {
	  sprintf(lbl, "%s %d", "SIA", i);
	  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	  sprintf(lbl, "%s %d", "CME", i);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_axpbycz(&CME, &SIA, &RIA, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
	  dpd_file2_close(&CME);
	  dpd_file2_close(&SIA);
	  sprintf(lbl, "%s %d", "Sia", i);
	  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
	  sprintf(lbl, "%s %d", "Cme", i);
	  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	  dpd_file2_axpbycz(&Cme, &Sia, &Ria, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
	  dpd_file2_close(&Cme);
	  dpd_file2_close(&Sia);
	}
      }

      if(params.ref == 0) norm = sqrt(2.0 * dpd_file2_dot_self(&RIA));
      else norm = sqrt(dpd_file2_dot_self(&RIA) + dpd_file2_dot_self(&Ria));

      if (pf) fprintf(outfile,"%22d%15.10lf%11.2e%12.2e",k+1,lambda[k],
		      lambda[k]-lambda_old[k], norm); 

      if ( (norm > eom_params.residual_tol_SS) ||
	   (fabs(lambda[k]-lambda_old[k]) > eom_params.eval_tol) ) {
	if (pf) fprintf(outfile,"%7s\n","N");

	if(params.ref == 0) precondition_SS_RHF(&RIA, lambda[k]);
	else precondition_SS(&RIA, &Ria, lambda[k]);

	if(params.ref == 0) {
	  norm = sqrt(2.0 * dpd_file2_dot_self(&RIA));
	  dpd_file2_scm(&RIA, 1.0/norm);
	}
	else {
	  norm = sqrt(dpd_file2_dot_self(&RIA) + dpd_file2_dot_self(&Ria));
	  dpd_file2_scm(&RIA, 1.0/norm);
	  dpd_file2_scm(&Ria, 1.0/norm);
	}

	if(params.ref == 0) schmidt_add_SS_RHF(&RIA, &numCs);
	else schmidt_add_SS(&RIA, &Ria, &numCs);
      }
      else {
	if (pf) fprintf(outfile,"%7s\n","Y");
	++num_converged;
	converged[k] = 1;
      }

      fflush(outfile);
    }

    dpd_file2_close(&RIA);
    if(params.ref) dpd_file2_close(&Ria);

    for (i=0;i<num_roots;++i) lambda_old[i] = lambda[i];
    free(lambda);

    if ( (iter == 2) || (L > eom_params.vectors_per_root_SS * eom_params.rpi[C_irr])) {
      restart_SS(alpha, L, eom_params.rpi[C_irr], C_irr);
      L = num_roots = eom_params.rpi[C_irr];
      keep_going = 1;
      free_block(alpha);
    }
    else {
      if (numCs > L) {
	keep_going = 1;
	free_block(alpha);
	L = numCs;
      }
    }
    ++iter;
  }

  if (pf) fprintf(outfile,"\nLowest eigenvalues of HBar Singles-Singles Block\n");
  if (pf) fprintf(outfile,"Root      Excitation Energy         Total Energy\n");
  if (pf) fprintf(outfile,"           (eV)     (cm^-1)            (au)\n");
  if (pf) for (i=0;i<num_roots;++i)
    if (pf)   fprintf(outfile,"%4d%12.3lf%12.2lf%20.10lf\n",i+1,
		      lambda_old[i]* _hartree2ev, lambda_old[i]* _hartree2wavenumbers,
		      lambda_old[i]+moinfo.eref+moinfo.ecc);

  free(lambda_old);
  /* collapse solutions to one vector each */
  restart_SS(alpha, L, eom_params.rpi[C_irr], C_irr);
  free_block(alpha);
  if (pf) fprintf(outfile,"\n");
  fflush(outfile);

  return;
}

void precondition_SS(dpdfile2 *RIA, dpdfile2 *Ria, double eval)
{
  dpdfile2 DIA, Dia;
  int h, nirreps, i, j, a, b, ij, ab;
  double tval;
  int C_irr = 0; /* to be determined by RIA irrep */

  nirreps = RIA->params->nirreps;

  dpd_file2_mat_init(RIA);
  dpd_file2_mat_rd(RIA);
  dpd_file2_init(&DIA, EOM_D, C_irr, 0, 1, "DIA");
  dpd_file2_mat_init(&DIA);
  dpd_file2_mat_rd(&DIA);
  for(h=0; h < nirreps; h++)
    for(i=0; i < RIA->params->rowtot[h]; i++)
      for(a=0; a < RIA->params->coltot[h]; a++) {
	tval = eval - DIA.matrix[h][i][a];
	if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
      }
  dpd_file2_mat_wrt(RIA);
  dpd_file2_mat_close(RIA);
  dpd_file2_close(&DIA);

  dpd_file2_mat_init(Ria);
  dpd_file2_mat_rd(Ria);
  dpd_file2_init(&Dia, EOM_D, C_irr, 0, 1, "Dia");
  dpd_file2_mat_init(&Dia);
  dpd_file2_mat_rd(&Dia);
  for(h=0; h < nirreps; h++)
    for(i=0; i < Ria->params->rowtot[h]; i++)
      for(a=0; a < Ria->params->coltot[h]; a++) {
	tval = eval - Dia.matrix[h][i][a];
	if (fabs(tval) > 0.0001) Ria->matrix[h][i][a] /= tval;
      }
  dpd_file2_mat_wrt(Ria);
  dpd_file2_mat_close(Ria);
  dpd_file2_close(&Dia);

  return;
}

void precondition_SS_RHF(dpdfile2 *RIA, double eval)
{
  dpdfile2 DIA;
  int h, nirreps, i, j, a, b, ij, ab;
  double tval;
  int C_irr = 0; /* to be determined by RIA irrep */

  nirreps = RIA->params->nirreps;

  dpd_file2_mat_init(RIA);
  dpd_file2_mat_rd(RIA);
  dpd_file2_init(&DIA, EOM_D, C_irr, 0, 1, "DIA");
  dpd_file2_mat_init(&DIA);
  dpd_file2_mat_rd(&DIA);
  for(h=0; h < nirreps; h++)
    for(i=0; i < RIA->params->rowtot[h]; i++)
      for(a=0; a < RIA->params->coltot[h]; a++) {
	tval = eval - DIA.matrix[h][i][a];
	if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
      }
  dpd_file2_mat_wrt(RIA);
  dpd_file2_mat_close(RIA);
  dpd_file2_close(&DIA);

  return;
}

void schmidt_add_SS(dpdfile2 *RIA, dpdfile2 *Ria, int *numCs)
{
  double dotval;
  double norm;
  int i, I;
  dpdfile2 Cme, CME;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  dpdbuf4 CMnEf_buf;
  char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];
  int C_irr = 0; /* to be determined by RIA irrep */

  for (i=0; i<*numCs; i++) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dotval  = dpd_file2_dot(RIA, &CME);

    sprintf(Cme_lbl, "%s %d", "Cme", i);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dotval += dpd_file2_dot(Ria, &Cme);

    dpd_file2_axpy(&CME, RIA, -1.0*dotval, 0);
    dpd_file2_axpy(&Cme, Ria, -1.0*dotval, 0);
    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }

  norm = dpd_file2_dot_self(RIA);
  norm += dpd_file2_dot_self(Ria);

  if (norm < eom_params.schmidt_add_residual_tol) {
    return;
  }
  else {
    scm_C1(RIA, Ria, 1.0/norm);
    /*
    fprintf(outfile, "Norm of residual RIA before scaling: %20.10f\n", norm);
    fprintf(outfile, "Norm of residual RIA to be added: %20.10f\n", dpd_file2_dot_self(RIA));
    fprintf(outfile, "Norm of residual Ria to be added: %20.10f\n", dpd_file2_dot_self(Ria));
    */
    sprintf(CME_lbl, "%s %d", "CME", *numCs);
    sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
    dpd_file2_copy(RIA, EOM_CME, CME_lbl);
    dpd_file2_copy(Ria, EOM_Cme, Cme_lbl);
    ++(*numCs);
  }

  return;
}

void schmidt_add_SS_RHF(dpdfile2 *RIA, int *numCs)
{
  double dotval;
  double norm;
  int i, I;
  dpdfile2 CME;
  char CME_lbl[32], Cme_lbl[32];
  int C_irr = 0; /* to be determined by RIA irrep */

  for (i=0; i<*numCs; i++) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dotval  = 2.0 * dpd_file2_dot(RIA, &CME);
    dpd_file2_axpy(&CME, RIA, -1.0*dotval, 0);
    dpd_file2_close(&CME);
  }

  norm = sqrt(2.0 * dpd_file2_dot_self(RIA));

  if (norm < eom_params.schmidt_add_residual_tol) {
    return;
  }
  else {
    dpd_file2_scm(RIA, 1.0/norm);
/*
    fprintf(outfile, "Norm of residual RIA before scaling: %20.10f\n", norm);
    fprintf(outfile, "Norm of residual RIA to be added: %20.10f\n", 2.0 * dpd_file2_dot_self(RIA));
*/
    sprintf(CME_lbl, "%s %d", "CME", *numCs);
    dpd_file2_copy(RIA, EOM_CME, CME_lbl);
    sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
    dpd_file2_copy(RIA, EOM_Cme, Cme_lbl);  /* to be removed after sigma spin-adaptation */
    ++(*numCs);
  }

  return;
}

void restart_SS(double **alpha, int L, int num, int C_irr) {
  int i,j,I;
  char lbl[20];
  dpdfile2 C1, CME, Cme, CME2, Cme2;
  dpdbuf4 C2, CMNEF, Cmnef, CMnEf;
  double norm, dotval;

  for (I=1;I<num;++I) {
    for (i=0; i<I; i++) {
      dotval = 0.0;
      for (j=0;j<L;++j) {
	dotval += alpha[j][i] * alpha[j][I];
      }
      for (j=0; j<L; j++) alpha[j][I] -= dotval * alpha[j][i];
    }
    dotval = 0.0;
    for (j=0;j<L;++j) dotval += alpha[j][I] * alpha[j][I];
    norm = sqrt(dotval);
    for (j=0;j<L;++j) alpha[j][I] = alpha[j][I]/norm;
  }


  for (i=0; i<num; ++i) {
    sprintf(lbl, "%s %d", "CME", L+i);
    dpd_file2_init(&C1, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_scm(&C1, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "CME", j);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_axpy(&CME, &C1, alpha[j][i], 0);
      dpd_file2_close(&CME);
    }
    dpd_file2_close(&C1);

    sprintf(lbl, "%s %d", "Cme", L+i);
    dpd_file2_init(&C1, EOM_Cme, C_irr, 0, 1, lbl);
    dpd_file2_scm(&C1, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "Cme", j);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
      dpd_file2_axpy(&Cme, &C1, alpha[j][i],0);
      dpd_file2_close(&Cme);
    }
    dpd_file2_close(&C1);
  }

  for (i=0; i<num; ++i) {
    sprintf(lbl, "%s %d", "CME", L+i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_copy(&CME, EOM_CME, lbl);
    dpd_file2_close(&CME);
    sprintf(lbl, "%s %d", "Cme", L+i);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    dpd_file2_copy(&Cme, EOM_Cme, lbl);
    dpd_file2_close(&Cme);

  }
  return;
}

