/*
**
** DETCAS
**
** Orbital rotation program for determinant configuration interaction
** wavefunctions evaluated using the DETCI program
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
**
** Modified 05/28/99 by CDS to move BLAS interface to libqt
**
*/

#include "globaldefs.h"
#include <math.h>

/* C INCLUDE FILES */

extern "C" {
   #include <stdlib.h>
   #include <stdio.h>
   #include <string.h>
   #include <libipv1/ip_lib.h>
   #include <libqt/qt.h>
   #include <libciomr/libciomr.h>
   #include <libfile30/file30.h>
   #include "globals.h"
   #include "setup_io.h"

   extern void get_mo_info(void);
   extern void get_parameters(void);
   extern void print_parameters(void);
   extern void read_integrals(void);
   extern void read_density_matrices(void);
   extern void read_lagrangian(void);
   extern void form_independent_pairs(void);
   extern void read_thetas(int npairs);
   extern void write_thetas(int npairs);
   extern int  read_ref_orbs(void);
   extern int  write_ref_orbs(void);
   extern void read_cur_orbs(void);
   extern void form_F_act(void);
   extern int  diis(int veclen, double *vec, double *errvec);
   extern void get_mat_block(double **src, double **dst, int dst_dim,
                             int dst_offset, int *dst2src);
   extern void calc_dE_dT(int n, double **dEU, int npairs, int *ppair, 
                          int *qpair, double *theta, double *dET);
   extern void form_diag_mo_hess(int npairs, int *ppair, int *qpair, 
                                 double *F_core, double *tei, double **opdm, 
                                 double *tpdm, double *F_act, int firstact, 
                                 int lastact, double *hess);
   extern void form_diag_mo_hess2(int npairs, int *ppair, int *qpair, 
                                 double *F_core, double *tei, double **opdm, 
                                 double *tpdm, double *F_act, int firstact, 
                                 int lastact, double *hess);
   extern void calc_orb_step(int npairs, double *grad, double *hess_diag, 
                             double *theta);
   extern int print_step(int npairs, int steptype);
   extern void postmult_by_U(int irrep, int dim, double **mo_coeffs,
                             int npairs, int *p_arr, int *q_arr, 
                             double *theta_arr);
   extern void premult_by_U(int irrep, int dim, double **mo_coeffs,
                            int npairs, int *p_arr, int *q_arr, 
                            double *theta_arr);
   extern void cleanup(void);
}


/* C++ INCLUDE FILES */

#include <iostream.h>
#include "indpairs.h"


/* C++ FUNCTION PROTOTYPES */

void parse_cmdline(int argc, char *argv[]);
void title(void);
void quote(void);
void init_ioff(void);
void calc_gradient(void);
void calc_hessian(void);
void scale_gradient(void);
int check_conv(void);
int take_step(void);
void rotate_orbs(void);


/* GLOBAL VARIABLES (other modules load these via globals.h) */
extern "C" {
  struct calcinfo CalcInfo;
  struct params Params;
  int *ioff;
  FILE *infile, *outfile;
}
IndepPairs IndPairs;



/* MAIN ROUTINE */

int main(int argc, char *argv[])
{
  int converged = 0;
  int num_pairs = 0;
  int steptype = 0;

  Params.print_lvl = 1;

  parse_cmdline(argc, argv);   /* check for any command-line arguments     */
  init_io();                   /* open input and output files              */
  get_parameters();            /* get running params (convergence, etc)    */
  init_ioff();                 /* set up the ioff array                    */
  title();                     /* print program identification             */

  if (Params.print_lvl) print_parameters();

  get_mo_info();               /* read DOCC, SOCC, frozen, nbfso, etc      */
  read_integrals();            /* get the 1 and 2 elec MO integrals        */
  read_density_matrices();
  read_lagrangian();

  form_independent_pairs();
  num_pairs = IndPairs.get_num_pairs();

  read_thetas(num_pairs);
  if (Params.print_lvl > 2)
    IndPairs.print_vec(CalcInfo.theta_cur, "\n\tRotation Angles:", outfile);

  if (!read_ref_orbs()) {
    read_cur_orbs();
    write_ref_orbs();
    zero_arr(CalcInfo.theta_cur, num_pairs);
    write_thetas(num_pairs);
  }

  form_F_act();
  calc_hessian();
  calc_gradient();
  converged = check_conv();
  scale_gradient();

  if (!converged) {
    steptype = take_step();
    rotate_orbs();
    write_thetas(num_pairs);
  }
  else
    steptype = 0;

  print_step(num_pairs, steptype);

  if (Params.print_lvl) quote();
  cleanup();
  close_io();
  return(converged);
}


/*
** parse_cmdline(): Figure out what the command-line arguments are.
**
*/
void parse_cmdline(int argc, char *argv[])
{
  int i;

  strcpy(Params.ofname, "output.dat");
  for (i=1; i<argc; i++) {
    if (strcmp(argv[i], "-quiet") == 0) Params.print_lvl = 0;
    else if (strcmp(argv[i], "-o") == 0) {
      if (strlen(argv[i+1]) < PARM_OUTFILE_MAX)
        strcpy(Params.ofname, argv[i+1]);
      else {
        fprintf(stderr, "detcas: output file name too long!\n");
      }
    }
  }

}


/*
** init_ioff(): Set up the ioff array for quick indexing
*/
void init_ioff(void)
{
  int i;

  /* set offsets for ij-type canonical ordering */
  ioff = (int *) malloc (IOFF_MAX * sizeof(int)) ;
  ioff[0] = 0;
  for (i = 1; i < IOFF_MAX ; i++) ioff[i] = ioff[i-1] + i;
}



/*
** title(): Function prints a program identification
*/
void title(void)
{
  if (Params.print_lvl) {
   fprintf(outfile,"\n");
   fprintf(outfile,"*******************************************************\n");
   fprintf(outfile,"                      D E T C A S \n");
   fprintf(outfile,"\n");
   fprintf(outfile,"                   C. David Sherrill\n") ;
   fprintf(outfile,"                     April 27 1998\n") ;
   fprintf(outfile,"*******************************************************\n");
   fprintf(outfile,"\n\n\n");
   }
  else {
   fprintf(outfile, 
     "\nD E T C A S: C. David Sherrill, April 27 1998\n");
   }
  fflush(outfile);
}


void quote(void)
{
  fprintf(outfile,"\n\t\t \"Good bug ... dead bug\" \n\n");
  fprintf(outfile,"\t\t\t - Ed Valeev\n\n");
  fflush(outfile);
}



void form_independent_pairs(void)
{

  IndPairs.set(CalcInfo.nirreps, MAX_RAS_SPACES, CalcInfo.ras_opi,
               CalcInfo.ras_orbs, CalcInfo.fzc_orbs, CalcInfo.fzv_orbs, 
               CalcInfo.frozen_docc, CalcInfo.frozen_uocc, CalcInfo.ci2relpitz, 
               Params.fci);

  if (Params.print_lvl > 3) IndPairs.print(outfile);

}


/*
** calc_gradient()
**
** This function calculates the MO gradient from the MO Lagrangian
**
*/
void calc_gradient(void)
{
  int pair, npair, h, ir_npairs, ir_norbs, offset;
  double *ir_mo_grad, **ir_lag, *ir_theta_cur, value, rms;
  int *parr, *qarr, *ir_ppair, *ir_qpair;

  npair = IndPairs.get_num_pairs();
  parr  = IndPairs.get_p_ptr();
  qarr  = IndPairs.get_q_ptr();

  CalcInfo.mo_grad = init_array(npair);

  /*
  calc_grad_1(npair, parr, qarr, CalcInfo.lag, CalcInfo.mo_grad);
  calc_grad_2(npair, parr, qarr, CalcInfo.onel_ints, CalcInfo.twoel_ints, 
              CalcInfo.opdm, CalcInfo.tpdm, CalcInfo.F_act, 
              (CalcInfo.num_cor_orbs + CalcInfo.num_fzc_orbs), 
              CalcInfo.npop, CalcInfo.mo_grad); 
  */

  // scratch array for dEdTheta, big enough for any irrep
  ir_mo_grad = init_array(npair);
  
  // calculate dEdU, then dEdTheta
  for (h=0,offset=0; h<CalcInfo.nirreps; h++) {

    // Setup for this irrep
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    ir_norbs = CalcInfo.orbs_per_irr[h];
    if (h>0) offset += CalcInfo.orbs_per_irr[h-1];
    if (!ir_npairs) continue;
    ir_ppair = IndPairs.get_ir_prel_ptr(h);
    ir_qpair = IndPairs.get_ir_qrel_ptr(h);
    ir_lag = block_matrix(ir_norbs, ir_norbs);
    get_mat_block(CalcInfo.lag, ir_lag, ir_norbs, offset, CalcInfo.pitz2ci);

    if (Params.print_lvl > 3) {
      fprintf(outfile, "Irrep %d of lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, outfile);
    }

    ir_theta_cur = IndPairs.get_irrep_vec(h, CalcInfo.theta_cur); 

    // Need to mult the Lagrangian by 2 to get dEdU
    C_DSCAL(ir_norbs*ir_norbs, 2.0, ir_lag[0], 1);

    if (Params.print_lvl > 3) {
      fprintf(outfile, "Irrep %d of 2 * lagrangian:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, outfile);
    }

    // Calc dEdU
    premult_by_U(h, CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
                 ir_ppair, ir_qpair, ir_theta_cur);

    if (Params.print_lvl > 3) {
      fprintf(outfile, "dE/dU:\n", h);
      print_mat(ir_lag, ir_norbs, ir_norbs, outfile);
    }

    // Calculate dEdTheta
    calc_dE_dT(CalcInfo.orbs_per_irr[h], ir_lag, ir_npairs,
               ir_ppair, ir_qpair, ir_theta_cur, ir_mo_grad);

    // Put dEdTheta into the large gradient array
    IndPairs.put_irrep_vec(h, ir_mo_grad, CalcInfo.mo_grad);
    delete [] ir_theta_cur;
    free_block(ir_lag); 
  }


  rms = 0.0;
  for (pair=0; pair<npair; pair++) {
    value =  CalcInfo.mo_grad[pair];
    rms += value * value;
  }

  if (Params.print_lvl > 2) 
    IndPairs.print_vec(CalcInfo.mo_grad, "\n\tOrbital Gradient:", outfile);

  rms = sqrt(rms);
  CalcInfo.mo_grad_rms = rms;

  if (Params.print_lvl)
    fprintf(outfile, "\n\tRMS Orbital Gradient: %6.4E\n", rms);
}


/*
** calc_hessian()
**
** This function calculates an approximate MO Hessian from the 
** Fock matrix intermediates and two-electron integrals.
**
** C. David Sherrill
** April 1998
*/
void calc_hessian(void)
{

  int npairs, *ppair, *qpair, ncore;

  npairs = IndPairs.get_num_pairs();
  ppair = IndPairs.get_p_ptr();
  qpair = IndPairs.get_q_ptr();

  /* allocate space for diagonal elements of MO Hessian */
  CalcInfo.mo_hess_diag = init_array(npairs);

  /* Now calculate the approximate diagonal MO Hessian */
  ncore = CalcInfo.num_fzc_orbs + CalcInfo.num_cor_orbs;

  if (strcmp(Params.hessian, "DIAG") == 0) {
    form_diag_mo_hess2(npairs, ppair, qpair, CalcInfo.onel_ints, 
                      CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm, 
                      CalcInfo.F_act, ncore, CalcInfo.npop, 
                      CalcInfo.mo_hess_diag);
    }
  else if (strcmp(Params.hessian, "APPROX_DIAG") == 0) {
    form_diag_mo_hess(npairs, ppair, qpair, CalcInfo.onel_ints, 
                      CalcInfo.twoel_ints, CalcInfo.opdm, CalcInfo.tpdm, 
                      CalcInfo.F_act, ncore, CalcInfo.npop, 
                      CalcInfo.mo_hess_diag);
    }
  else {
    fprintf(outfile, "(detcas): Unrecognized Hessian option %s\n", Params.hessian);
    }
 

  if (Params.print_lvl > 3)
    IndPairs.print_vec(CalcInfo.mo_hess_diag, 
                       "\n\tApproximate Diagonal MO Hessian:", outfile);

}



/*
** scale_gradient()
**
** Scales the orbital gradient by the approximate orbital Hessian
*/
void scale_gradient(void)
{
  int pair, npairs;
  double rms, value;

  npairs = IndPairs.get_num_pairs();

  CalcInfo.theta_step = init_array(npairs);

  // All this actually does is scale the gradient by the Hessian
  if (Params.scale_grad)
    calc_orb_step(npairs, CalcInfo.mo_grad, CalcInfo.mo_hess_diag,
                  CalcInfo.theta_step);
  else {
    for (pair=0; pair<npairs; pair++) 
      CalcInfo.theta_step[pair] = CalcInfo.mo_grad[pair];
  }
    

  if (Params.print_lvl > 3)  
    IndPairs.print_vec(CalcInfo.theta_step,"\n\tScaled Orbital Grad:", outfile);

  rms = 0.0;
  for (pair=0; pair<npairs; pair++) {
    value =  CalcInfo.theta_step[pair];
    rms += value * value;
  }

  rms = sqrt(rms);
  CalcInfo.scaled_mo_grad_rms = rms;
 
  if (Params.print_lvl)
    fprintf(outfile, "\n\tScaled RMS Orbital Gradient: %6.4E\n", rms);

  if (Params.scale_step != 1.0) {
    for (pair=0; pair<npairs; pair++) 
      CalcInfo.theta_step[pair] *= Params.scale_step;
  }

}


/*
** take_step()
**
** This function takes a step in orbital rotation (theta) space
**
** Returns: type of step taken; 1=regular (Newton-Raphson), 2=diis
*/
int take_step(void)
{
  int npairs, pair, took_diis;
  
  npairs = IndPairs.get_num_pairs();

  for (pair=0; pair<npairs; pair++) 
    CalcInfo.theta_cur[pair] += CalcInfo.theta_step[pair];

  if (CalcInfo.iter >= Params.diis_start) 
    took_diis = diis(npairs, CalcInfo.theta_cur, CalcInfo.theta_step);
    //took_diis = diis(npairs, CalcInfo.theta_cur, CalcInfo.mo_grad);
  else 
    took_diis = 0;

  if (!took_diis) {
    if (Params.print_lvl) 
      fprintf(outfile, "Taking regular step\n");
  }

  return(took_diis+1);

}


/*
** rotate_orbs()
**
** Rotate the orbitals, irrep by irrep
*/
void rotate_orbs(void)
{
  double *ir_theta;
  int h, pair, ir_norbs, ir_npairs, *ir_ppair, *ir_qpair;

  // First, we need to come up with Theta vectors for each irrep
  ir_theta = init_array(IndPairs.get_num_pairs());  // always big enough
  for (h=0; h<CalcInfo.nirreps; h++) {
    ir_npairs = IndPairs.get_ir_num_pairs(h);
    if (ir_npairs) {
      ir_norbs = CalcInfo.orbs_per_irr[h];
      ir_theta = IndPairs.get_irrep_vec(h, CalcInfo.theta_cur);
      ir_ppair  = IndPairs.get_ir_prel_ptr(h);
      ir_qpair  = IndPairs.get_ir_qrel_ptr(h);
      
      if (Params.print_lvl > 3) {
        fprintf(outfile, "Thetas for irrep %d\n", h);
        for (pair=0; pair<ir_npairs; pair++) {
          fprintf(outfile, "Pair (%2d,%2d) = %12.6lf\n",
                  ir_ppair[pair], ir_qpair[pair], ir_theta[pair]);
        }
        fprintf(outfile, "\n");
        fflush(outfile);
      }

      /* print new coefficients */
      if (Params.print_mos) {
        fprintf(outfile, "\n\tOld molecular orbitals for irrep %s\n", 
          CalcInfo.labels[h]);
        print_mat(CalcInfo.mo_coeffs[h], ir_norbs, ir_norbs, outfile);
      }

      postmult_by_U(h, ir_norbs, CalcInfo.mo_coeffs[h], ir_npairs, 
                    ir_ppair, ir_qpair, ir_theta);

      /* print new coefficients */
      if (Params.print_mos) {
        fprintf(outfile, "\n\tNew molecular orbitals for irrep %s\n", 
          CalcInfo.labels[h]);
        print_mat(CalcInfo.mo_coeffs[h], ir_norbs, ir_norbs, outfile);
      }


      /* write the new block of MO coefficients to file30 */
      file30_init();
      file30_wt_blk_scf(CalcInfo.mo_coeffs[h], h);
      file30_close();

      delete [] ir_theta;
    }
  }


}


/*
** check_conv
**
** Check the summary file to see if we've converged
**
** Returns: 1 if converged, otherwise 0
*/
int check_conv(void)
{
  FILE *sumfile;
  char sumfile_name[] = "file14.dat";
  char comment[MAX_COMMENT];
  int i, entries, iter, nind;
  double rmsgrad, scaled_rmsgrad, energy, energy_last;
  int converged_energy=0, converged_grad=0, last_converged=0;
  double conv_rms_grad, conv_e;

  sumfile = fopen(sumfile_name, "r");

  if (sumfile == NULL) {
    CalcInfo.iter = 0;
    return(0);
  }

  if (fscanf(sumfile, "%d", &entries) != 1) {
    fprintf(outfile,"(print_step): Trouble reading num entries in file %s\n",
            sumfile_name);
    fclose(sumfile);
    CalcInfo.iter = 0;
    return(0);
  } 

  CalcInfo.iter = entries;
  for (i=0; i<entries; i++) {
    fscanf(sumfile, "%d %d %lf %lf %lf %s", &iter, &nind, &scaled_rmsgrad,
           &rmsgrad, &energy_last, comment);
  }
  fclose(sumfile);

  file30_init();
  energy = file30_rd_ecorr();
  file30_close();

  /* check for convergence */
  conv_rms_grad = pow(10.0, -(Params.rms_grad_convergence));
  conv_e = pow(10.0, -(Params.energy_convergence));
  if (rmsgrad < conv_rms_grad) converged_grad = 1;
  if (fabs(energy_last - energy) < conv_e)
    converged_energy = 1;
  if (strstr(comment, "CONV") != NULL)
    last_converged = 1;

  if (converged_grad && converged_energy && !last_converged) {
    fprintf(outfile, "\n\t*** Calculation Converged ***\n");
    return(1);
  }
  else {
    fprintf(outfile, "\n\t... calculation continuing ...\n");
    return(0);
  }

}

