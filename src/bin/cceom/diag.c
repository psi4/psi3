/*
 *   diag() diagonalized Hbar with Davidson-Liu algorithm to obtain
 *    right-hand eigenvector and eigenvalue
 */

#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"
#include <physconst.h>

extern void test_dpd();
void rzero(int C_irr);
void rzero_rhf(int C_irr);
void init_S1(int index, int irrep);
void init_S2(int index, int irrep);
void init_C1(int i, int C_irr);
void init_C2(int index, int irrep);
extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);
extern void scm_C1(dpdfile2 *CME, dpdfile2 *Cme, double a);
extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
    dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);
extern void scm_C2(dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);
extern void restart(double **alpha, int L, int num, int irrep, int ortho);
extern void precondition(dpdfile2 *RIA, dpdfile2 *Ria,
    dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, double eval);
extern void precondition_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, double eval);
void form_diagonal(int irrep);
extern void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria,
    dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep);
extern void schmidt_add_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int *numCs, int irrep);
void c_clean(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
void sigmaSS(int index, int irrep);
void sigmaSD(int index, int irrep);
void sigmaDS(int index, int irrep);
void sigmaDD(int index, int irrep);
void diagSS(int irrep);
void hbar_extra(void);
void hbar_norms(void);
extern void sort_C(int index, int irrep);

double local_G1_dot(dpdfile2 *, dpdfile2 *);
double local_G2_dot(dpdbuf4 *, dpdbuf4 *);
void local_filter_T1_nodenom(dpdfile2 *);
void local_filter_T2_nodenom(dpdbuf4 *);

void read_guess(int);

void diag(void) {
  dpdfile2 CME, Cme, SIA, Sia, RIA, Ria, DIA, Dia, tIA, tia, LIA, Lia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, SIJAB, Sijab, SIjAb, RIJAB, Rijab, RIjAb, RIjbA;
  dpdbuf4 CMnEf1, CMnfE1, CMnfE, CMneF, C2;
  char lbl[32];
  int num_converged, *converged, keep_going, already_sigma;
  int irrep, numCs, iter, lwork, info;
  int get_right_ev = 1, get_left_ev = 0;
  int L,h,i,j,k,a,nirreps,errcod,C_irr;
  double norm, tval, **G, *work, *evals_complex, **alpha, **evectors_left;
  double *lambda, *lambda_old;
  int num_vecs;
  double ra, rb, r2aa, r2bb, r2ab;

  hbar_extra(); /* sort hbar matrix elements for sigma equations */
#ifdef EOM_DEBUG
  /* hbar_norms(); */
#endif

  fprintf(outfile,"Symmetry of ground state: %s\n", moinfo.labels[moinfo.sym]);
  /* loop over symmetry of C's */
  for (C_irr=0; C_irr<moinfo.nirreps; ++C_irr) {
    already_sigma = 0;
    iter = 0;
    keep_going = 1;
    num_converged = 0;

    if (eom_params.cs_per_irrep[C_irr] == 0) continue;
    fprintf(outfile,"Symmetry of excited state: %s\n", moinfo.labels[moinfo.sym ^ C_irr]);
    fprintf(outfile,"Symmetry of right eigenvector: %s\n",moinfo.labels[C_irr]);

    /* Store approximate diagonal elements of Hbar */
    form_diagonal(C_irr);

    if(params.local) {
      local_guess();
      if(local.do_singles) diagSS(C_irr);
    }
    else {
      if(!strcmp(eom_params.guess,"SINGLES")) {
	/* Diagonalize Hbar-SS to obtain initial CME and Cme guess */
	fprintf(outfile,"Obtaining initial guess from singles-singles block of Hbar...");
	diagSS(C_irr);
	if (!eom_params.print_singles) fprintf(outfile,"Done.\n\n");
      }
      else if(!strcmp(eom_params.guess,"INPUT")) {
	read_guess(C_irr);
      }
    }

    /* Uncomment to print out initial guesses */
    /*
      for (i=0;i<eom_params.cs_per_irrep[C_irr]; ++i) {
      sprintf(lbl, "%s %d", "CME", i);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_print(&CME,outfile);
      dpd_file2_close(&CME);
      if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Cme", i);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
      dpd_file2_print(&Cme,outfile);
      dpd_file2_close(&Cme);
      }
      }
    */

    /* Setup and zero initial C2 and S2 vector to go with Hbar_SS */
    for (i=0;i<eom_params.cs_per_irrep[C_irr];++i) {
      /* init_S1(i, C_irr); gets done at first iteration anyway */
      init_C2(i, C_irr);
      /* init_S2(i, C_irr); */
    }

#ifdef EOM_DEBUG
    check_sum("reset", 0, 0); /* reset checksum */ 
#endif

    converged = init_int_array(eom_params.cs_per_irrep[C_irr]);
    lambda_old = init_array(eom_params.cs_per_irrep[C_irr]);
    L = eom_params.cs_per_irrep[C_irr];

    while ((keep_going == 1) && (iter < eom_params.max_iter)) {
      fprintf(outfile,"Iter=%-4d L=%-4d", iter+1, L); fflush(outfile);
      keep_going = 0;
      numCs = L;
      num_converged = 0;

      for (i=already_sigma;i<L;++i) {
        /* Form a zeroed S vector for each C vector
	   SIA and Sia do get overwritten by sigmaSS
	   so this may only be necessary for debugging */
        init_S1(i, C_irr);
        init_S2(i, C_irr);

        sort_C(i, C_irr);

        /* Computing sigma vectors */
#ifdef TIME_CCEOM
        timer_on("sigma all");
        timer_on("sigmaSS"); sigmaSS(i,C_irr); timer_off("sigmaSS");
        timer_on("sigmaSD"); sigmaSD(i,C_irr); timer_off("sigmaSD");
        timer_on("sigmaDS"); sigmaDS(i,C_irr); timer_off("sigmaDS");
        timer_on("sigmaDD"); sigmaDD(i,C_irr); timer_off("sigmaDD");
        timer_off("sigma all");
#else
        sigmaSS(i,C_irr);
        sigmaSD(i,C_irr);
        sigmaDS(i,C_irr);
        sigmaDD(i,C_irr);
#endif

        /* Cleaning out sigma vectors for open-shell cases  */
        if (params.eom_ref > 0) {
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
          check_sum("reset",0,0);
          sprintf(lbl, "Total sigma%d norm", i);
          check_sum(lbl, i, C_irr);
#endif
          c_clean(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);

#ifdef EOM_DEBUG
          check_sum("reset",0,0);
          sprintf(lbl, "Total sigma%d norm af clean", i);
          check_sum(lbl, i, C_irr);
#endif

          dpd_file2_close(&SIA);
          dpd_file2_close(&Sia);
          dpd_buf4_close(&SIJAB);
          dpd_buf4_close(&Sijab);
          dpd_buf4_close(&SIjAb);
          fflush(outfile);
        }
      }

      already_sigma = L;

      /* Form G = C'*S matrix */
      G = block_matrix(L,L);
      for (i=0;i<L;++i) {

        if(params.eom_ref == 0) {
          /* Spin-adapt C */
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);

          sprintf(lbl, "%s %d", "CMnEf", i);
          dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
          dpd_buf4_copy(&CMnEf, EOM_TMP, "CMnEf");
          dpd_buf4_sort(&CMnEf, EOM_TMP, pqsr, 0, 5, "CMnfE");
          dpd_buf4_close(&CMnEf);

          dpd_buf4_init(&CMnEf1, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnEf");
          dpd_buf4_scm(&CMnEf1, 2.0);
          dpd_buf4_init(&CMnfE1, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
          dpd_buf4_axpy(&CMnfE1, &CMnEf1, -1.0);
          dpd_buf4_close(&CMnfE1);
          dpd_buf4_close(&CMnEf1);

          /* dpd_file2_init(&CME, EOM_TMP, C_irr, 0, 1, "CME");*/
          dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnEf");
        }
        else {

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

        }

        for (j=0;j<L;++j) {

          if(params.eom_ref == 0) {
            sprintf(lbl, "%s %d", "SIA", j);
            dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
	    tval = 2.0 * dpd_file2_dot(&CME, &SIA);
            sprintf(lbl, "%s %d", "SIjAb", j);
            dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
	    tval += dpd_buf4_dot(&CMnEf, &SIjAb);
            dpd_file2_close(&SIA);
            dpd_buf4_close(&SIjAb);
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
          }
          G[i][j] = tval;
        }

        dpd_file2_close(&CME);
        dpd_buf4_close(&CMnEf);
        if (params.eom_ref > 0) {
          dpd_file2_close(&Cme);
          dpd_buf4_close(&CMNEF);
          dpd_buf4_close(&Cmnef);
        }
      } /* end build of G */
      /* mat_print(G, L, L, outfile); */
      /* Diagonalize G Matrix */ 
      lambda = init_array(L);        /* holds real part of eigenvalues of G */
      alpha = block_matrix(L,L);     /* will hold eigenvectors of G */
      dgeev_eom(L, G, lambda, alpha);

      eigsort(lambda, alpha, L);

      /* eivout(alpha, lambda, L, L, outfile);*/

      free_block(G);

      /* Compute Residual vectors */
      dpd_file2_init(&RIA, EOM_R, C_irr, 0, 1, "RIA");
      dpd_buf4_init(&RIjAb, EOM_R, C_irr, 0, 5, 0, 5, 0, "RIjAb");
      if (params.eom_ref > 0) {
        dpd_file2_init(&Ria, EOM_R, C_irr, 0, 1, "Ria");
        dpd_buf4_init(&RIJAB, EOM_R, C_irr, 2, 7, 2, 7, 0, "RIJAB");
        dpd_buf4_init(&Rijab, EOM_R, C_irr, 2, 7, 2, 7, 0, "Rijab");
      }
      fprintf(outfile,"  Root    EOM Energy     Delta E   Res. Norm    Conv?\n");
      for (k=0;k<eom_params.cs_per_irrep[C_irr];++k) {
        /* rezero residual vector for each root */
        dpd_file2_scm(&RIA, 0.0);
        dpd_buf4_scm(&RIjAb, 0.0);
        if (params.eom_ref > 0) {
          dpd_file2_scm(&Ria, 0.0);
          dpd_buf4_scm(&RIJAB, 0.0);
          dpd_buf4_scm(&Rijab, 0.0);
        }

        converged[k] = 0;
        for (i=0;i<L;++i) {
          sprintf(lbl, "%s %d", "SIA", i);
          dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
          sprintf(lbl, "%s %d", "CME", i);
          dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
          dpd_file2_axpbycz(&CME, &SIA, &RIA, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_file2_close(&CME);
          dpd_file2_close(&SIA);

          sprintf(lbl, "%s %d", "CMnEf", i);
          dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
          sprintf(lbl, "%s %d", "SIjAb", i);
          dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
          dpd_buf4_axpbycz(&CMnEf, &SIjAb, &RIjAb, -1.0*lambda[k]*alpha[i][k], alpha[i][k], 1.0);
          dpd_buf4_close(&CMnEf);
          dpd_buf4_close(&SIjAb);
          if (params.eom_ref > 0) {
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
          }
        }

	if(params.eom_ref == 0) precondition_RHF(&RIA, &RIjAb, lambda[k]);
	else precondition(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, lambda[k]);

        if (params.eom_ref == 0) {
          dpd_buf4_sort(&RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
          dpd_buf4_init(&RIjbA, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");
          norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
        }
        else norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);

#ifdef EOM_DEBUG
        fprintf(outfile,"Norm of residual vector %d  bf clean %18.13lf\n",k,norm);
#endif
        /* necessary? c_clean(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb); 
	   norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb); */
#ifdef EOM_DEBUG
        fprintf(outfile,"Norm of residual vector %d af clean %18.13lf\n",k,norm);
#endif

        fprintf(outfile,"%22d%15.10lf%11.2e%12.2e",k+1,lambda[k],
		lambda[k]-lambda_old[k], norm);

        /* Check for convergence and add new vector if not converged */
        if ( (norm > eom_params.residual_tol) ||
	     (fabs(lambda[k]-lambda_old[k]) > eom_params.eval_tol) ) {
          fprintf(outfile,"%7s\n","N");

	    /*  if(params.eom_ref == 0) precondition_RHF(&RIA, &RIjAb, lambda[k]);
	        else precondition(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, lambda[k]); */

        if(params.eom_ref == 0) {

	    /* if(params.local) {
	      local_filter_T1(&RIA, 0);
	      local_filter_T2(&RIjAb, 0);
	      } */

            dpd_buf4_sort(&RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
            dpd_buf4_init(&RIjbA, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "RIjbA");

            norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
            dpd_buf4_close(&RIjbA);

            dpd_file2_scm(&RIA, 1.0/norm);
            dpd_buf4_scm(&RIjAb, 1.0/norm);
          }
          else {
            norm = norm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb);
            scm_C(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, 1.0/norm);
          }

#ifdef EOM_DEBUG
          fprintf(outfile,"Norm of residual vector af preconditioning %18.13lf\n",norm);
#endif

          if(params.eom_ref == 0) schmidt_add_RHF(&RIA, &RIjAb, &numCs, C_irr);
          else schmidt_add(&RIA, &Ria, &RIJAB, &Rijab, &RIjAb, &numCs, C_irr);
        }
        else {
          fprintf(outfile,"%7s\n","Y");
          ++num_converged;
          converged[k] = 1;
        }
      }

      dpd_file2_close(&RIA);
      dpd_buf4_close(&RIjAb);
      if (params.eom_ref > 0) {
        dpd_file2_close(&Ria);
        dpd_buf4_close(&RIJAB);
        dpd_buf4_close(&Rijab);
      }

      for (i=0;i<eom_params.cs_per_irrep[C_irr];++i) lambda_old[i] = lambda[i];
      free(lambda);

      /* restart with new B vectors if there are too many */
      if (L > eom_params.vectors_per_root * eom_params.cs_per_irrep[C_irr]) {
        restart(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr, 1);
        L = eom_params.cs_per_irrep[C_irr];
        keep_going = 1;
        /* already_sigma = 0; */
        already_sigma = L;
      }
      else {
        /* If any new vectors were added, then continue */
        if (numCs > L) {
          keep_going = 1;
          L = numCs;
        }
      }

      ++iter;
      /* If we're all done then collapse to one vector per root */
      if ( (keep_going == 0) && (iter < eom_params.max_iter) ) {
        fprintf(outfile,"Collapsing to only %d vectors.\n", eom_params.cs_per_irrep[C_irr]);
        restart(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr, 0);
      }
      free_block(alpha);
    }

    if (C_irr == eom_params.prop_sym^moinfo.sym) {
      /* Copy desired root to CC_RAMPS file */
      fprintf(outfile,"Copying C for root %d to CC_RAMPS\n",eom_params.prop_root);
      sprintf(lbl, "%s %d", "CME", eom_params.prop_root-1);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_copy(&CME, CC_RAMPS, "RIA");
      dpd_file2_close(&CME);
      sprintf(lbl, "%s %d", "CMnEf", eom_params.prop_root-1);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_copy(&CMnEf, CC_RAMPS, "RIjAb");
      dpd_buf4_close(&CMnEf);
  
      if(params.eom_ref > 0) {
        sprintf(lbl, "%s %d", "Cme", eom_params.prop_root-1);
        dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
        dpd_file2_copy(&Cme, CC_RAMPS, "Ria");
        dpd_file2_close(&Cme);
        sprintf(lbl, "%s %d", "CMNEF", eom_params.prop_root-1);
        dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
        dpd_buf4_copy(&CMNEF, CC_RAMPS, "RIJAB");
        dpd_buf4_close(&CMNEF);
        sprintf(lbl, "%s %d", "Cmnef", eom_params.prop_root-1);
        dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
        dpd_buf4_copy(&Cmnef, CC_RAMPS, "Rijab");
        dpd_buf4_close(&Cmnef);
      }
    }

    /* Print out results */
    fprintf(outfile,"\nProcedure converged for %d roots.\n",num_converged);

    if (num_converged == eom_params.cs_per_irrep[C_irr]) { }
    else if (iter == eom_params.max_iter) {
      fprintf(outfile,"\nMaximum number of iterations exceeded,");
      fprintf(outfile,"so not all roots converged.\n");
    }
    else {
      fprintf(outfile,"\nAlgorithm failure: No vectors could be added,");
      fprintf(outfile,"though not all roots converged.\n");
    }

    if (num_converged > 0) {
      fprintf(outfile,"\nFinal Energetic Summary for Converged Roots of Irrep %s\n",
	      moinfo.labels[moinfo.sym^C_irr]);
      fprintf(outfile," Root           Excitation Energy           Total Energy\n");
      fprintf(outfile,"          (eV)     (cm^-1)    (au)            (au)\n");
      for (i=0;i<eom_params.cs_per_irrep[C_irr];++i) {
        if (converged[i] == 1) {
          fprintf(outfile,"%4d %10.3lf %10.1lf %14.10lf  %15.10lf\n",i+1,
		  lambda_old[i]* _hartree2ev, lambda_old[i]* _hartree2wavenumbers,
		  lambda_old[i], lambda_old[i]+moinfo.eref+moinfo.ecc);

	  if(params.eom_ref == 0) {

	    /* Print out largest few components of this root */
	    fprintf(outfile, "\nLargest components of excited wave function #%d:\n", i);

	    sprintf(lbl, "%s %d", "CME", i);
	    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
	    amp_write_T1(&CME, 5, outfile);
	    dpd_file2_close(&CME);

	    sprintf(lbl, "%s %d", "CMnEf", i);
	    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
	    amp_write_T2(&CMnEf, 5, outfile);
	    dpd_buf4_close(&CMnEf);
	    fprintf(outfile, "\n");
	  }

        }
      }
      psio_write_entry(CC_INFO, "CCEOM Energy",
		       (char *) &(lambda_old[eom_params.prop_root-1]), sizeof(double));

      i = moinfo.sym ^ C_irr;
      psio_write_entry(CC_INFO, "CCEOM State Irrep", (char *) &i, sizeof(int));

      /*      fprintf(outfile,"\nCCEOM energy %.10lf and state irrep %d written to CC_INFO.\n",
	      lambda_old[eom_params.prop_root-1], i); */
    }
    fprintf(outfile,"\n");

    /* remove all temporary files */
  free(lambda_old);
  free(converged);
    for(i=CC_TMP; i<CC_RAMPS; i++) psio_close(i,0);
    for(i=CC_TMP; i<CC_RAMPS; i++) psio_open(i,0);
  }

  if (params.eom_ref == 0) rzero_rhf(eom_params.prop_sym^moinfo.sym);
  else rzero(eom_params.prop_sym^moinfo.sym);

  return;
}

/* zeroes ith CME (and Cme) vectors on disk */
void init_C1(int i, int C_irr ){
  dpdfile2 CME, Cme;
  char lbl[32];
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_scm(&CME, 0.0);
    dpd_file2_close(&CME);
  }
  else {
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    if (params.eom_ref == 1) dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
    else if (params.eom_ref == 2) dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
    scm_C1(&CME, &Cme, 0.0);
    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
}

/* zeroes ith SIA (and Sia) vectors on disk */
void init_S1(int i, int C_irr) {
  dpdfile2 SIA, Sia;
  char lbl[32];
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    dpd_file2_scm(&SIA, 0.0);
    dpd_file2_close(&SIA);
  }
  else {
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    if (params.eom_ref == 1) dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
    else if (params.eom_ref == 2) dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
    scm_C1(&SIA, &Sia, 0.0);
    dpd_file2_close(&SIA);
    dpd_file2_close(&Sia);
  }
}

/* zeroes ith CMnEf (+ CMNEF + Cmnef) on disk */
void init_C2(int i, int C_irr) {
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&CMnEf, 0.0);
    dpd_buf4_close(&CMnEf);
  }
  else {
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);

    sprintf(lbl, "%s %d", "Cmnef", i);
    if (params.eom_ref == 1)
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);

    sprintf(lbl, "%s %d", "CMnEf", i);
    if (params.eom_ref == 1)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);

    /* scm_C2(&CMNEF, &Cmnef, &CMnEf, 0.0); */
    dpd_buf4_scm(&CMNEF, 0.0);
    dpd_buf4_scm(&Cmnef, 0.0);
    dpd_buf4_scm(&CMnEf, 0.0);

    dpd_buf4_close(&CMNEF);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_close(&CMnEf);
  }
}


void init_S2(int i, int C_irr) {
  dpdbuf4 SIJAB, Sijab, SIjAb;
  char lbl[32];
  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&SIjAb, 0.0);
    dpd_buf4_close(&SIjAb);
  }
  else {
    sprintf(lbl, "%s %d", "SIJAB", i);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Sijab", i);
    if (params.eom_ref == 1) 
      dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
    else if (params.eom_ref == 2) 
      dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, lbl);

    sprintf(lbl, "%s %d", "SIjAb", i);
    if (params.eom_ref == 1)
      dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);
    scm_C2(&SIJAB, &Sijab, &SIjAb, 0.0);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&SIjAb);
  }
}

