/*
*   diag() diagonalized Hbar with Davidson-Liu algorithm to obtain
*    right-hand eigenvector and eigenvalue
*/

#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"
#include <physconst.h>

void init_S1(int index, int irrep);
void init_S2(int index, int irrep);
void init_C2(int index, int irrep);
double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
  dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
void restart(double **alpha, int L, int num, int irrep, int ortho);
void precondition(dpdfile2 *RIA, dpdfile2 *Ria,
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, double eval);
void form_diagonal(int irrep);
void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria,
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep);
void c_clean(dpdfile2 *CME, dpdfile2 *Cme,
  dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
void sigmaSS(int index, int irrep);
void sigmaSD(int index, int irrep);
void sigmaDS(int index, int irrep);
void sigmaDD(int index, int irrep);
void diagSS(int irrep);
void hbar_extra(void);

void diag(void) {
  dpdfile2 Fmi, FMI, Fae, FAE, Fme, FME;
  dpdfile2 CME, Cme, SIA, Sia, RIA, Ria;
  dpdbuf4 CMNEF, Cmnef, CMnEf, SIJAB, Sijab, SIjAb, RIJAB, Rijab, RIjAb;
  dpdfile2 tIA, tia;
  dpdbuf4 tIJAB, tijab, tIjAb, W;
  dpdfile2 DIA, Dia;
  dpdbuf4 DIJAB, Dijab, DIjAb;

  char lbl[32], lbl2[32];
  int num_converged = 0, *converged, keep_going = 1, already_sigma = 0;
  int irrep, numCs, iter = 0, lwork, info;
  int get_right_ev = 1, get_left_ev = 0;
  int L,h,i,j,k,a,nirreps,errcod,C_irr;
  double norm, tval, **G, *work, *evals_complex, **alpha, **evectors_left;
  double *lambda, *lambda_old;

  hbar_extra();

  for (irrep=0; irrep<moinfo.nirreps; ++irrep) {
    if (eom_params.rpi[irrep] == 0) continue; /* no roots of this irrep desired */
  
    for (C_irr=0;C_irr<moinfo.nirreps;++C_irr) /* determine symmetry of C */
      if ( (moinfo.sym ^ C_irr) == irrep ) break;

    fprintf(outfile,"\nSymmetry of right eigenvector: %s\n",moinfo.labels[C_irr]);

    if (C_irr != 0) {
      fprintf(outfile,"\nUnable to handle asymmetric R amplitudes for now\n");
      continue;
    }

    /* Store approximate diagonal elements of Hbar */
    form_diagonal(C_irr);

    /* Diagonalize Hbar-SS to obtain initial CME and Cme guess */
    diagSS(C_irr);

    /* Setup initial C2 and S2 vector to go with Hbar_SS */
    for (i=0;i<eom_params.rpi[C_irr];++i) {
      init_C2(i, C_irr);
      init_S2(i, C_irr);
    }
    check_sum("sending zero", 0, 0); /* reset checksum by sending zero vector */

    converged = init_int_array(eom_params.rpi[irrep]);
    lambda_old = init_array(eom_params.rpi[irrep]);
    L = eom_params.rpi[C_irr];

    while ((keep_going == 1) && (iter < eom_params.max_iter)) {
      fprintf(outfile,"Iter=%-4d L=%-4d", iter+1, L); fflush(outfile);
      keep_going = 0;
      numCs = L;
      num_converged = 0;

      for (i=already_sigma;i<L;++i) {
	   /* Form a zeroed S vector for each C vector */
	   /* SIA and Sia do get overwritten by sigmaSS */
           /* so this may only be necessary for debugging */
        init_S1(i, C_irr);
        init_S2(i, C_irr);

	/* Make copies of current C vectors */
	/* Copy used in WmbejDD */
	sprintf(lbl, "%s %d", "CMNEF", i);
	dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
	dpd_buf4_sort(&CMNEF, EOM_TMP, prqs, 10, 10, "CMENF");
	dpd_buf4_close(&CMNEF);
	sprintf(lbl, "%s %d", "Cmnef", i);
	dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 0, 5, 2, 7, 0, lbl);
	dpd_buf4_sort(&Cmnef, EOM_TMP, prqs, 10, 10, "Cmenf");
	dpd_buf4_close(&Cmnef);
	sprintf(lbl, "%s %d", "CMnEf", i);
	dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
	dpd_buf4_sort(&CMnEf, EOM_TMP, prqs, 10, 10, "CMEnf");
	/* Copy used in WmnieSD */
	dpd_buf4_sort(&CMnEf, EOM_TMP, qprs, 0, 5, "CmNEf");
	/* Copy of current C vector used in WmnieSD and WabefDD */
	dpd_buf4_sort(&CMnEf, EOM_TMP, pqsr, 0, 5, "CMneF");
	dpd_buf4_close(&CMnEf);
	/* Copy used in WmbejDD */
	dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMneF");
	dpd_buf4_sort(&CMnEf, EOM_TMP, prqs, 10, 10, "CMenF");
	dpd_buf4_close(&CMnEf);
	dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNEf");
	dpd_buf4_sort(&CMnEf, EOM_TMP, prqs, 10, 10, "CmENf");
	dpd_buf4_close(&CMnEf);
	/* Copy used in FDD, FSD, WamefSD, WmnefDD, WmnieSD */
	dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNEf");
	dpd_buf4_sort(&CMnEf, EOM_TMP, pqsr, 0, 5, "CmNeF");
	dpd_buf4_close(&CMnEf);
	/* Copy used in WmbejDD */

       dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNeF");
       dpd_buf4_sort(&CMnEf, EOM_TMP, prqs, 10, 10, "CmeNF");
       dpd_buf4_close(&CMnEf);

/*
	dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMEnf");
	dpd_buf4_sort(&CMnEf, EOM_TMP, rspq, 10, 10, "CmeNF");
	dpd_buf4_close(&CMnEf);
*/

	/* Computing sigma vectors */
	sigmaSS(i,C_irr);
/*
	sigmaSD(i,C_irr);
	sigmaDS(i,C_irr);
	sigmaDD(i,C_irr);
*/

	/* Cleaning out sigma vectors */
	sprintf(lbl, "%s %d", "SIA", i);
	dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	sprintf(lbl, "%s %d", "Sia", i);
	dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
	sprintf(lbl, "%s %d", "SIJAB", i);
	dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
	sprintf(lbl, "%s %d", "Sijab", i);
	dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
	sprintf(lbl, "%s %d", "SIjAb", i);
	dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

#ifdef EOM_DEBUG
	norm = norm_C(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);
	fprintf(outfile,"Norm of sigma %d bf clean: %18.13lf\n",i,norm);
#endif
/*
	c_clean(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);
*/

#ifdef EOM_DEBUG
	norm = norm_C(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);
	fprintf(outfile,"Norm of sigma %d af clean: %18.13lf\n",i,norm);
#endif

	dpd_file2_close(&SIA);
	dpd_file2_close(&Sia);
	dpd_buf4_close(&SIJAB);
	dpd_buf4_close(&Sijab);
	dpd_buf4_close(&SIjAb);
        fflush(outfile);
      }

      already_sigma = L;

      /* Form G = C'*S matrix */
      G = block_matrix(L,L);
      for (i=0;i<L;++i) {
	sprintf(lbl, "%s %d", "CME", i);
	dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	sprintf(lbl, "%s %d", "Cme", i);
	dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	sprintf(lbl, "%s %d", "CMNEF", i);
	dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
	sprintf(lbl, "%s %d", "Cmnef", i);
	dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
	sprintf(lbl, "%s %d", "CMnEf", i);
	dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
	for (j=0;j<L;++j) {
	  sprintf(lbl, "%s %d", "SIA", j);
	  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	  tval = dpd_file2_dot(&CME, &SIA);
	  dpd_file2_close(&SIA);
	  sprintf(lbl, "%s %d", "Sia", j);
	  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
	  tval += dpd_file2_dot(&Cme, &Sia);
	  dpd_file2_close(&Sia);
	  sprintf(lbl, "%s %d", "SIJAB", j);
	  dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
	  tval += dpd_buf4_dot(&CMNEF, &SIJAB);
	  dpd_buf4_close(&SIJAB);
	  sprintf(lbl, "%s %d", "Sijab", j);
	  dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
	  tval += dpd_buf4_dot(&Cmnef, &Sijab);
	  dpd_buf4_close(&Sijab);
	  sprintf(lbl, "%s %d", "SIjAb", j);
	  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
	  tval += dpd_buf4_dot(&CMnEf, &SIjAb);
	  dpd_buf4_close(&SIjAb);
	  G[i][j] = tval;
	}
	dpd_file2_close(&Cme);
	dpd_file2_close(&CME);
	dpd_buf4_close(&CMNEF);
	dpd_buf4_close(&Cmnef);
	dpd_buf4_close(&CMnEf);
      }


      /* Diagonalize G Matrix */

      lambda = init_array(L);        /* holds real part of eigenvalues of G */
      alpha = block_matrix(L,L);     /* will hold eigenvectors of G */

      dgeev_eom(L, G, lambda, alpha);

/*
      dgeev_call(&get_left_ev, &get_right_ev, &L, &(G[0][0]), &L, &(lambda[0]),
		 &(evals_complex[0]), &(evectors_left[0][0]), &L, &(alpha[0][0]),
		 &L, &(work[0]), &lwork, &info);
*/

      eigsort(lambda, alpha, L);

      free_block(G);

      /* Compute Residual vectors */
      dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
      dpd_file2_init(&Ria, EOM_R, C_irr, 0, 1, "Ria");
      dpd_buf4_init(&RIJAB, EOM_R, C_irr, 2, 7, 2, 7, 0, "RIJAB");
      dpd_buf4_init(&Rijab, EOM_R, C_irr, 2, 7, 2, 7, 0, "Rijab");
      dpd_buf4_init(&RIjAb, EOM_R, C_irr, 0, 5, 0, 5, 0, "RIjAb");
      fprintf(outfile,"  Root   EOM Energy     Delta E    Res. Norm    Conv?\n");
      for (k=0;k<eom_params.rpi[irrep];++k) {
        /* rezero residual vector for each root */
        scm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, 0.0);
	converged[k] = 0;
	for (i=0;i<L;++i) {
          sprintf(lbl, "%s %d", "SIA", i);
	  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "CME", i);
	  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	  dpd_file2_axpbycz(&CME, &SIA, &RIA,
            -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_file2_close(&CME);
          dpd_file2_close(&SIA);

          sprintf(lbl, "%s %d", "Sia", i);
	  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "Cme", i);
	  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
	  dpd_file2_axpbycz(&Cme, &Sia, &Ria,
            -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_file2_close(&Cme);
          dpd_file2_close(&Sia);

          sprintf(lbl, "%s %d", "CMNEF", i);
          dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "SIJAB", i);
          dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
	  dpd_buf4_axpbycz(&CMNEF, &SIJAB, &RIJAB, 
            -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_buf4_close(&CMNEF);
          dpd_buf4_close(&SIJAB);

          sprintf(lbl, "%s %d", "Cmnef", i);
          dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
          sprintf(lbl, "%s %d", "Sijab", i);
          dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
	  dpd_buf4_axpbycz(&Cmnef, &Sijab, &Rijab,
            -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_buf4_close(&Cmnef);
          dpd_buf4_close(&Sijab);

          sprintf(lbl, "%s %d", "CMnEf", i);
          dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
          sprintf(lbl, "%s %d", "SIjAb", i);
          dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_axpbycz(&CMnEf, &SIjAb, &RIjAb,
            -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_buf4_close(&CMnEf);
          dpd_buf4_close(&SIjAb);
	}

	norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);

#ifdef EOM_DEBUG
	fprintf(outfile,"Norm of residual vector bf clean %18.13lf\n",norm);
#endif
	/* necessary? */
/*
	c_clean(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);
*/
	norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);
#ifdef EOM_DEBUG
	fprintf(outfile,"Norm of residual vector af clean %18.13lf\n",norm);
#endif

	fprintf(outfile,"%21d%15.10lf%12.2e%12.2e",k+1,lambda[k],
		lambda[k]-lambda_old[k], norm);

	/* Check for convergence and add new vector if not converged */
	if ( (norm > eom_params.residual_tol) ||
	     (fabs(lambda[k]-lambda_old[k]) > eom_params.eval_tol) ) {
	  fprintf(outfile,"%7s\n","N");
	  precondition(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, lambda[k]);
	  norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);
	  scm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, 1.0/norm);
	  schmidt_add(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, &numCs, irrep);
	}
	else {
          fprintf(outfile,"%7s\n","Y");
          ++num_converged;
          converged[k] = 1;
	}
      }
      dpd_file2_close(&RIA);
      dpd_file2_close(&Ria);
      dpd_buf4_close(&RIJAB);
      dpd_buf4_close(&Rijab);
      dpd_buf4_close(&RIjAb);

      for (i=0;i<eom_params.rpi[irrep];++i) lambda_old[i] = lambda[i];
      free(lambda);

      /* restart with new B vectors if there are too many */
      if (L > eom_params.vectors_per_root * eom_params.rpi[irrep]) {
        restart(alpha, L, eom_params.rpi[irrep], C_irr, 1);
        L = eom_params.rpi[irrep];
        keep_going = 1;
        already_sigma = 0;
      }
      else {
        /* If any new vectors were added, then continue */
        if (numCs > L) {
	  keep_going = 1;
	  L = numCs;
        }
      }

      ++iter;
      /* If we're all done then collapse vectors (before alpha is freed) */
      if ( (keep_going == 0) && (iter < eom_params.max_iter) ) {
        fprintf(outfile,"We're gonna quit so calling restart.\n");
        restart(alpha, L, eom_params.rpi[irrep], C_irr, 0);
      }
      free_block(alpha);
    }

    /* Copy desired root to CC_RAMPS file */
    fprintf(outfile,"Copying C %d to CC_RAMPS\n",eom_params.prop_root-1);
    sprintf(lbl, "%s %d", "CME", eom_params.prop_root-1);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, CC_OEI, "RIA");
    dpd_file2_close(&CME);

    sprintf(lbl, "%s %d", "Cme", eom_params.prop_root-1);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
    dpd_file2_copy(&Cme, CC_OEI, "Ria");
    dpd_file2_close(&Cme);

    sprintf(lbl, "%s %d", "CMNEF", eom_params.prop_root-1);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&CMNEF, CC_RAMPS, "RIJAB");
    dpd_buf4_close(&CMNEF);

    sprintf(lbl, "%s %d", "Cmnef", eom_params.prop_root-1);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&Cmnef, CC_RAMPS, "Rijab");
    dpd_buf4_close(&Cmnef);

    sprintf(lbl, "%s %d", "CMnEf", eom_params.prop_root-1);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_copy(&CMnEf, CC_RAMPS, "RIjAb");
    dpd_buf4_close(&CMnEf);

    /* Print out results */
    fprintf(outfile,"\nProcedure converged for %d roots.\n",num_converged);

    if (num_converged == eom_params.rpi[irrep]) { }
    else if (iter == eom_params.max_iter) {
      fprintf(outfile,"\nMaximum number of iterations exceeded,");
      fprintf(outfile,"so not all roots converged.\n");
    }
    else {
      fprintf(outfile,"\nAlgorithm failure: No vectors could be added,");
      fprintf(outfile,"though not all roots converged.\n");
    }

    if (num_converged > 0) {
      fprintf(outfile,"\nFinal Energetic Summary for Converged Roots of Irrep %d\n", C_irr);
      fprintf(outfile,"Root      Excitation Energy         Total Energy\n");
      fprintf(outfile,"           (eV)      (cm^-1)            (au)\n");
      for (i=0;i<eom_params.rpi[irrep];++i) {
        if (converged[i] == 1) {
	  fprintf(outfile,"%4d%12.5lf%12.1lf%20.10lf\n",i+1,
		  lambda_old[i]* _hartree2ev, lambda_old[i]* _hartree2wavenumbers,
		  lambda_old[i]+moinfo.eref+moinfo.ecc);
        }
      }
      psio_write_entry(CC_INFO, "CCEOM Energy", (char *) &(lambda_old[eom_params.prop_root-1]),
		       sizeof(double));
      fprintf(outfile,"\nThe CCEOM energy %.10lf is written to CC_INFO.\n",
	      lambda_old[eom_params.prop_root-1]);

      psio_write_entry(CC_INFO, "CCEOM Irrep", (char *) &irrep, sizeof(int));

      fprintf(outfile,"\nThe CCEOM irrep %d is written to CC_INFO.\n");

    }
    fprintf(outfile,"\n");
    free(lambda_old);
  }
  return;
}


void restart(double **alpha, int L, int num, int C_irr, int ortho) {
  int i,I,j,h;
  char lbl[20];
  dpdfile2 C1, CME, Cme;
  dpdbuf4 C2, CMNEF, Cmnef, CMnEf;
  double dotval, norm;

  /* Orthonormalize alpha[1] through alpha[num] */
if (ortho) {
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
}

  /* Form restart vectors Ci = Sum_j(alpha[j][i]*Cj) */
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
       dpd_file2_axpy(&Cme, &C1, alpha[j][i], 0);
       dpd_file2_close(&Cme);
    }
    dpd_file2_close(&C1);
    sprintf(lbl, "%s %d", "CMNEF", L+i);
    dpd_buf4_init(&C2, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
       sprintf(lbl, "%s %d", "CMNEF", j);
       dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
       dpd_buf4_axpy(&CMNEF, &C2, alpha[j][i]);
       dpd_buf4_close(&CMNEF);
    }
    dpd_buf4_close(&C2);
    sprintf(lbl, "%s %d", "Cmnef", L+i);
    dpd_buf4_init(&C2, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
       sprintf(lbl, "%s %d", "Cmnef", j);
       dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
       dpd_buf4_axpy(&Cmnef, &C2, alpha[j][i]);
       dpd_buf4_close(&Cmnef);
    }
    dpd_buf4_close(&C2);
    sprintf(lbl, "%s %d", "CMnEf", L+i);
    dpd_buf4_init(&C2, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
       sprintf(lbl, "%s %d", "CMnEf", j);
       dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
       dpd_buf4_axpy(&CMnEf, &C2, alpha[j][i]);
       dpd_buf4_close(&CMnEf);
    }
    dpd_buf4_close(&C2);
  }

  /* Copy restart vectors to beginning of file */
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
    sprintf(lbl, "%s %d", "CMNEF", L+i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_copy(&CMNEF, EOM_CMNEF, lbl);
    dpd_buf4_close(&CMNEF);
    sprintf(lbl, "%s %d", "Cmnef", L+i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_copy(&Cmnef, EOM_Cmnef, lbl);
    dpd_buf4_close(&Cmnef);
    sprintf(lbl, "%s %d", "CMnEf", L+i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_copy(&CMnEf, EOM_CMnEf, lbl);
    dpd_buf4_close(&CMnEf);
  }

  return;
}

void init_S1(int i, int C_irr) {
  int h;
  dpdfile2 SIA, Sia;
  char lbl[32];
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "Sia", i);
  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
  scm_C1(&SIA, &Sia, 0.0);
  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);
}

/*
void init_S1(int i, int C_irr) {
  int h;
  dpdfile2 SIA, Sia;
  char lbl[32];
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  dpd_file2_mat_init(&SIA);
  dpd_file2_mat_wrt(&SIA);
  dpd_file2_mat_close(&SIA);
  sprintf(lbl, "%s %d", "Sia", i);
  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
  dpd_file2_mat_init(&Sia);
  dpd_file2_mat_wrt(&Sia);
  dpd_file2_mat_close(&Sia);
  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);
}
*/

void init_S2(int i, int C_irr) {
  dpdbuf4 SIJAB, Sijab, SIjAb;
  char lbl[32];
  sprintf(lbl, "%s %d", "SIJAB", i);
  dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "Sijab", i);
  dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

  dpd_buf4_scm(&SIJAB, 0.0);
  dpd_buf4_scm(&Sijab, 0.0);
  dpd_buf4_scm(&SIjAb, 0.0);

  dpd_buf4_close(&SIJAB);
  dpd_buf4_close(&Sijab);
  dpd_buf4_close(&SIjAb);
}

/*
void init_S2(int i, int C_irr) {
  int h;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  char lbl[32];
  sprintf(lbl, "%s %d", "SIJAB", i);
  dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "Sijab", i);
  dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
  for(h=0; h < SIJAB.params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(&SIJAB, h);
    dpd_buf4_mat_irrep_wrt(&SIJAB, h);
    dpd_buf4_mat_irrep_close(&SIJAB, h);
  }
  for(h=0; h < Sijab.params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Sijab, h);
    dpd_buf4_mat_irrep_wrt(&Sijab, h);
    dpd_buf4_mat_irrep_close(&Sijab, h);
  }
  for(h=0; h < SIjAb.params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(&SIjAb, h);
    dpd_buf4_mat_irrep_wrt(&SIjAb, h);
    dpd_buf4_mat_irrep_close(&SIjAb, h);
  }
  dpd_buf4_close(&SIJAB);
  dpd_buf4_close(&Sijab);
  dpd_buf4_close(&SIjAb);
  return;
}
*/

void init_C2(int i, int C_irr) {
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];
  sprintf(lbl, "%s %d", "CMNEF", i);
  dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "Cmnef", i);
  dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "CMnEf", i);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_scm(&CMNEF, 0.0);
  dpd_buf4_scm(&Cmnef, 0.0);
  dpd_buf4_scm(&CMnEf, 0.0);
  dpd_buf4_close(&CMNEF);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_close(&CMnEf);
}


/*
void init_C2(int i, int C_irr) {
  int h;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];
  sprintf(lbl, "%s %d", "CMNEF", i);
  dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "Cmnef", i);
  dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "CMnEf", i);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
  for(h=0; h < CMNEF.params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(&CMNEF, h);
    dpd_buf4_mat_irrep_wrt(&CMNEF, h);
    dpd_buf4_mat_irrep_close(&CMNEF, h);
  }
  for(h=0; h < Cmnef.params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Cmnef, h);
    dpd_buf4_mat_irrep_wrt(&Cmnef, h);
    dpd_buf4_mat_irrep_close(&Cmnef, h);
  }
  for(h=0; h < CMnEf.params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(&CMnEf, h);
    dpd_buf4_mat_irrep_wrt(&CMnEf, h);
    dpd_buf4_mat_irrep_close(&CMnEf, h);
  }
  dpd_buf4_close(&CMNEF);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_close(&CMnEf);
}
*/
