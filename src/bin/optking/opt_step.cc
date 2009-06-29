/*! \file
    \ingroup OPTKING
    \brief OPT_STEP.CC takes geometry steps using
    gradients -- optking's default operation
*/

#include <cmath>
#include <cstdio>
#include <libchkpt/chkpt.h>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

namespace psi { namespace optking {

extern double *compute_q(internals &simples, salc_set &symm);
extern double **compute_B(internals &simples, salc_set &symm);
extern double **compute_G(double **B, int num_intcos, cartesians &carts);
extern int new_geom(cartesians &carts, internals &simples, salc_set &symm, double *dq,
    int print_to_geom_file, int restart_geom_file, 
    char *disp_label, int disp_num, int last_disp, double *return_geom);
void fconst_init(cartesians &carts, internals &simples, salc_set &symm);
extern void compute_zmat(cartesians &carts, int *unique_zvars);
extern void print_zmat(FILE *outfile, int *unique_zvars);
extern double **compute_H(internals &simples, salc_set &symm, double **P, cartesians &carts);
extern void opt_report(FILE *of);
extern void free_info(int nsimples);

inline double rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t)/(_hartree2J*1.0e18);
}
inline double nr_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(_hartree2J*1.0e18);
}

bool line_search(cartesians &carts, int num_ints, double *dq);

int opt_step(cartesians &carts, internals &simples, salc_set &symm) {

  int xyz, i,j,k,ii,a,b, dim, dim_carts, nsimples, constraint, cnt;
  double **B, **G, **G2, **G_inv, **H, **H_inv, **temp_mat, **u, **P;
  double DE_error, DE, DE_old, E_old, DE_new;
  double **C, **T, **T2, **T3, R_AB, theta_A, theta_B, tau, phi_A, phi_B;
  double *temp_arr2, *temp_arr, *masses, **geom, *forces, *coord, max_force, rms_force;
  double *f, *f_q, *dq, *dq_to_new_geom, *q, tval, tval2, tval3, tval4, scale, temp, *djunk;
  char *disp_label, *wfn, value_string[30], force_string[30];
  bool do_line_search = false, success;

  int rfo_root, rfo_num;
  double **rfo_mat, *lambda, *work, *rfo_old_evect, *rfo_u;
  double rfo_xnorm, rfo_g, rfo_h, rfo_t, rfo_eval;
  double nr_xnorm, nr_g, nr_h, *nr_u;

  // for analytic interfragment coordinates
  int isalc, I, simple_id, intco_type, sub_index, sub_index2, atom;
  double **geom_A, **geom_B, **weight_A, **weight_B;
  double inter_q[6];

  dim_carts = carts.get_natom()*3;
  dim = symm.get_num();
  djunk = new double[dim_carts];
  disp_label = new char[MAX_LINELENGTH];

  if (dim == 0) punt("No symmetric internal coordinates to optimize.\n");

  ip_string("WFN", &(wfn),0);
  fprintf(outfile,"\nCurrent %s energy before step   %20.12lf\n", wfn, carts.get_energy());
  free(wfn);

  open_PSIF();
  if (optinfo.iteration > 0) {
    sprintf(value_string,"Energy %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &E_old, sizeof(double));
    fprintf(outfile,"  Previous energy                %20.12lf\n", E_old);
    DE = carts.get_energy() - E_old;
    fprintf(outfile,"  Actual Delta(E) from last step %20.12lf\n", DE);
    sprintf(value_string,"DE prediction %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &(DE_old), sizeof(double));
    DE_error = (DE_old-DE)/fabs(DE_old);
    fprintf(outfile,"  Previous Predicted Delta(E)    %20.12lf  (%5.1lf%% error)\n", DE_old,DE_error*100);
  }
  close_PSIF();

  // searching for a minimum, supersized steps downward permitted
  if (!optinfo.ts && (DE_error < -1*optinfo.step_energy_limit) && fabs(DE_old) > 1.0e-9 ) {
    fprintf(outfile,"\n\tInsufficient energy drop observed.\n");
    do_line_search = true;
  }

  // searching for a TS, no supersized steps permitted
  if (optinfo.ts && (fabs(DE_error) > optinfo.step_energy_limit) && fabs(DE_old) > 1.0e-9 ) {
    fprintf(outfile,"\t\tEnergy deviated too much from projection.\n");
    do_line_search = true;
  }

  // *** Bad step - step back and do line search
  if (do_line_search) {
    dq = init_array(dim);
    if (! line_search(carts,dim,dq) ) {
      fprintf(outfile,"Unable to complete line search. Quitting.\n");
      exit_io();
      exit(PSI_RETURN_FAILURE);
    }
    f_q = init_array(dim);
    open_PSIF(); // read old forces from previous step - should be same
    sprintf(force_string,"Internal Forces %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, force_string, (char *) f_q, dim* sizeof(double));
    sprintf(value_string,"RMS force %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &rms_force, sizeof(double));
    sprintf(value_string,"MAX force %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) &max_force, sizeof(double));
    close_PSIF();
    // update internal coordinate data to new cartesan geometry
    coord = carts.get_coord();
    simples.compute_internals(carts.get_natom(), coord);
    simples.compute_s(carts.get_natom(), coord);
    q = compute_q(simples,symm); 
  }
  else {
  fprintf(outfile,"\nTaking geometry step number %d\n", optinfo.iteration+1);

  // compute forces in internal coordinates, f_q = G_inv B u f
  q = compute_q(simples,symm);
  B = compute_B(simples,symm); //print_mat(B, dim, dim_carts,outfile);
  G = compute_G(B,dim,carts);

  fprintf(outfile,"\nBuB^t ");
  G_inv = symm_matrix_invert(G,dim,1,optinfo.redundant);

  masses = carts.get_mass();
  u = mass_mat(masses);
  free(masses);

  f = carts.get_forces(); // in aJ/Ang //print_mat2(&f, 1, dim_carts, outfile);

  f_q = init_array(dim);
  temp_arr = init_array(dim);
  temp_arr2 = init_array(dim_carts);

  mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
  mmult(B,0,&temp_arr2,1,&temp_arr,1,dim,dim_carts,1,0);
  mmult(G_inv,0,&temp_arr,1,&f_q,1,dim,dim,1,0);

  // fprintf(outfile,"internal forces (G_inv B u f) in aJ/A\n");
  // print_mat2(&f_q, 1, dim, outfile);
  // for (i=0;i<dim;++i)
  //   f_q[i] = f_q[i] * sqrt(2) * _hartree2J * 1.0E18 / _bohr2angstroms;
  // fprintf(outfile, "%d %15.10lf\n", i, f_q[i] * _hartree2J / _bohr2angstroms * 1.0E18) ;

  // test by transforming f_q back to cartesian forces and compare
  // fprintf(outfile,"computed forces in cartesian coordinates aJ/Ang\n");
  // mmult(B,1,&f_q,1,&temp_arr2,1,dim_carts,dim,1,0);
  // print_mat2(&temp_arr2, 1, dim_carts, outfile);

  free_block(B);
  free(f);
  free(temp_arr);
  free(temp_arr2);
  free_block(u);

  // Setup projection matrix P = G * G- for redundant internals
  P = block_matrix(dim,dim);
  mmult(G,0,G_inv,0,P,0,dim,dim,dim,0); 
  free_block(G);
  free_block(G_inv);
  
  // Add constraints to projection matrix
    // C has 1's on diagonal for constraints and 0's elsewhere
  if (optinfo.constraints_present) {
    C = block_matrix(dim,dim);
    double i_overlap, j_overlap;
    for (b=0; b<optinfo.nconstraints; ++b) {
      constraint = simples.index_to_id(optinfo.constraints[b]);

      for (i=0; i<dim; ++i) {
        for (j=0; j<=i; ++j) {
          i_overlap = 0.0;
          for (k=0;k<symm.get_length(i);++k) {
            a = symm.get_simple(i,k);
            if ( a != constraint ) continue;
            i_overlap += symm.get_coeff(i,k) * symm.get_prefactor(i);
          }
          j_overlap = 0.0;
          for (k=0;k<symm.get_length(j);++k) {
            a = symm.get_simple(j,k);
            if ( a != constraint ) continue;
            j_overlap += symm.get_coeff(j,k) * symm.get_prefactor(j);
          }
          C[i][j] += ( i_overlap * j_overlap );
          C[j][i] = C[i][j];
        }
      }
    }

  // P = P' - P' C (CPC)^-1 C P'
    T = block_matrix(dim,dim);
    T2 = block_matrix(dim,dim);
    mmult(P,0,C,0,T,0,dim,dim,dim,0); 
    mmult(C,0,T,0,T2,0,dim,dim,dim,0); 
    T3 = symm_matrix_invert(T2, dim, 0, 1);

    mmult(C,0,P,0,T,0,dim,dim,dim,0); 
    mmult(T3,0,T,0,T2,0,dim,dim,dim,0); 
    mmult(C,0,T2,0,T3,0,dim,dim,dim,0); 
    mmult(P,0,T3,0,T2,0,dim,dim,dim,0); 
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j)
        P[i][j] -= T2[i][j];
    free_block(T);
    free_block(T2);
    free_block(T3);
    free_block(C);
  }

   // Project redundancies and constraints out of forces: f_q = P f_q'
   temp_arr = init_array(dim);
   mmult(P,0,&f_q,1,&temp_arr,1,dim, dim,1,0);
   for(i=0; i<dim; ++i)
     f_q[i] = temp_arr[i];
   free(temp_arr);

  // Compute RMS and MAX forces
  rms_force = 0.0; 
  max_force = fabs(f_q[0]);
  for (i=0;i<dim;++i) {
    rms_force += SQR(f_q[i]);
    if (fabs(f_q[i]) > max_force) max_force = fabs(f_q[i]);
  } 
  rms_force = sqrt(rms_force/((double) dim));

  // Write data to PSIF_OPTKING for later use
  open_PSIF();
  sprintf(value_string,"Internal Values %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) q, dim*sizeof(double));
  sprintf(force_string,"Internal Forces %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, force_string, (char *) f_q, dim* sizeof(double));
  coord = carts.get_coord();
  sprintf(value_string,"Cartesian Values %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) coord, dim_carts*sizeof(double));
  free(coord);
  tval = carts.get_energy();
  sprintf(value_string,"Energy %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &tval, sizeof(double));
  sprintf(value_string,"RMS force %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &rms_force, sizeof(double));
  sprintf(value_string,"MAX force %d", optinfo.iteration);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &max_force, sizeof(double));
  close_PSIF();

  printf("Energy: %15.10lf MAX force: %6.2e RMS force: %6.2e\n",carts.get_energy(),max_force, rms_force);
      
  if (max_force < optinfo.conv) {
    fprintf(outfile, "\n\t       *** Internal Coordinates and Forces *** \n");
    fprintf(outfile, "\t-----------------------------------------------------------\n");
    fprintf(outfile, "\t Coordinate    Value (Ang or Rad)   Force aJ/Ang or aJ/Rad \n");
    fprintf(outfile, "\t-----------------------------------------------------------\n");
    for (i=0;i<dim;++i)
      fprintf(outfile,"\t %6d %20.8lf %20.8lf \n", i+1,  q[i],  f_q[i] );
    fprintf(outfile, "\t-----------------------------------------------------------\n");
    fprintf(outfile,"\nMAX force is < %5.1e.  Optimization is complete.\n", optinfo.conv);
    ip_string("WFN", &(wfn),0);
    fprintf(outfile,"\nFinal %s energy is %15.10lf\n", wfn, carts.get_energy());
    free(wfn);

    optinfo.iteration += 1; // set for the NEXT step
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
    close_PSIF();

    opt_report(outfile);

    fprintf(stdout,"\n  OPTKING:  optimization is complete\n");
    fprintf(outfile,"The Optimized geometry in a.u.\n");
    carts.print(12,outfile,0,disp_label,0);
    fprintf(outfile,"\nThe Optimized geometry in Angstrom\n");
    carts.print(13,outfile,0,disp_label,0);
    fprintf(outfile, "\n");

    if (optinfo.zmat) {
      int *unique_zvars;
      unique_zvars = init_int_array(MAX_ZVARS);
      compute_zmat(carts, unique_zvars);
      print_zmat(outfile, unique_zvars);
      free(unique_zvars);
      fprintf(outfile,"\n");
    }
    free(f_q); free(q);
    fprintf(stdout,"\n  Returning code %d\n", PSI_RETURN_ENDLOOP);
    return(PSI_RETURN_ENDLOOP);
  } // end converged geometry

  fconst_init(carts, simples, symm); // makes sure some force constants are in PSIF_OPTKING

  H = compute_H(simples,symm,P,carts); // get updated and projected Hessian
  free_block(P);

  // *** standard Newton-Raphson search ***
  if (!optinfo.rfo) { // displacements from inverted, projected Hessian, H_inv f_q = dq
    dq = init_array(dim);
    H_inv = symm_matrix_invert(H,dim,0,0);
    mmult(H_inv,0,&f_q,1,&dq,1,dim,dim,1,0);
    free_block(H_inv);

    /* determine scale factor needed to keep step less than 10% of q if q big or less than 0.1
     if q not so big, hopefully better solution later */
    scale = 1.0;
    for (i=0;i<dim;++i) {
      if (fabs(q[i]) > optinfo.step_limit) { // intco is not very small
        if (fabs(dq[i]) > STEP_PERCENT*fabs(q[i])) { // step is too large
          temp = STEP_PERCENT*fabs(q[i])/fabs(dq[i]);
        }
        else
          temp = 1;
      }
      else { // intco is very small
        if (fabs(dq[i]) < optinfo.step_limit)
          temp = 1.0; // step is small enough
        else
          temp = optinfo.step_limit / fabs(dq[i]);
      }
      if (temp < scale){
        scale = temp;
      }
    } 
    fprintf(outfile,"\nScaling displacements by %lf\n",scale);
    for (i=0;i<dim;++i)
      dq[i] = dq[i] * scale;   

    // avoid piddly symmetry breaking
    for (i=0;i<dim;++i)
      if (fabs(dq[i]) < MIN_DQ_STEP) dq[i] = 0.0;

    // get norm |x| and unit vector in the step direction
    nr_u = init_array(dim);
    for (i=0; i<dim; ++i)
      nr_u[i] = dq[i];
    dot_arr(nr_u, nr_u, dim, &tval);
    nr_xnorm = sqrt(tval);
    normalize(&nr_u, 1, dim);
    
    // get gradient and hessian in step direction
    dot_arr(f_q, nr_u, dim, &nr_g);
    nr_g *= -1; // gradient not force
    nr_h = 0;
    for (i=0; i<dim; ++i) {
      dot_arr(H[i], nr_u, dim, &tval);
      nr_h += nr_u[i] * tval;
    }

    DE_new = nr_energy(nr_xnorm, nr_g, nr_h);
    fprintf(outfile,"Projected energy change by quadratic approximation: %20.10lf\n", DE_new);
    
    open_PSIF();
    sprintf(value_string,"DE prediction %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &DE_new, sizeof(double));
    sprintf(value_string,"Unit step U %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) nr_u, dim*sizeof(double));
    sprintf(value_string,"Xnorm %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &nr_xnorm, sizeof(double));
    sprintf(value_string,"Scalar gradient %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &nr_g, sizeof(double));
    sprintf(value_string,"Scalar hesian %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &nr_h, sizeof(double));
    close_PSIF();
  }
  else { // take RFO step

    // build and diagonalize RFO matrix
    rfo_mat = block_matrix(dim+1,dim+1);
    for (i=0; i<dim; ++i) {
      rfo_mat[dim][i] = - f_q[i];
      for (j=0; j<=i; ++j)
        rfo_mat[i][j] = H[i][j];
    }

    if(optinfo.print_hessian) {
      fprintf(outfile,"RFO mat\n");
      print_mat5(rfo_mat,dim+1,dim+1,outfile);
    }

    lambda = init_array(dim+1);
    work = init_array(3*(dim+1));
    C_DSYEV('V', 'U', dim+1 , rfo_mat[0], dim+1, lambda, work, 3*(dim+1));
    free(work);
    for (i=0; i<dim; ++i)
      lambda[i] /= (_hartree2J * 1.0e18); //aJ -> hartree

  fprintf(outfile,"RFO eigenvalues/lambdas\n");
  print_mat2(&lambda,1,dim+1,outfile);

  // I'm not sure about the intermediate normalization that is supposed to occur
  // scaling the evects to make the last element '1' and then "scaling of the eigenvector
  // by a factor alpha used to calculate mode displacements" (see page 54) amounts to
  // the same thing as just making sure the sign of the last entry of the evect is positive
  //  During the course of an optimization some evects may appear that are bogus leads (presumably)
  // due to a bad hessian - the root following can avoid them.  also some evects have a virtual
  // 0 in the last entry and so cannot be divided by it
  for (i=0; i<dim+1; ++i) {
    if (rfo_mat[i][dim] < 0)
      for (j=0;j<dim+1;++j) rfo_mat[i][j] *= -1;
  }

  //fprintf(outfile,"RFO evects\n");
  //print_mat5(rfo_mat,dim+1,dim+1,outfile);

    // *** choose which RFO eigenvector to use
    // if not root following, then use rfo_root'th lowest eigenvalue; default is 0 (lowest)
    if (!optinfo.rfo_follow_root) {
      rfo_root = optinfo.rfo_root;
    }
    else { // do root following
      open_PSIF();
      if (psio_tocscan(PSIF_OPTKING, "Previous RFO eigenvector") != NULL) {
        rfo_old_evect = init_array(dim+1);
        psio_read_entry(PSIF_OPTKING, "Previous RFO eigenvector",
          (char *) &(rfo_old_evect[0]), (dim+1)*sizeof(double));

        tval = 0;
        for (i=0; i<dim+1; ++i) {
          dot_arr(rfo_mat[i],rfo_old_evect,dim+1,&tval2);
          if (tval2 > tval) {
            tval = tval2;
            rfo_root = i;
           }
        }
        free(rfo_old_evect);
        fprintf(outfile,"RFO vector %d has maximal overlap with previous step\n",rfo_root+1);
        psio_write_entry(PSIF_OPTKING, "Previous RFO eigenvector",
          (char *) rfo_mat[rfo_root], (dim+1)*sizeof(double));
      }
      else {
        psio_write_entry(PSIF_OPTKING, "Previous RFO eigenvector",
          (char *) &(rfo_mat[optinfo.rfo_root][0]), (dim+1)*sizeof(double));
        fprintf(outfile,"Using initial RFO vector %d to follow.\n",optinfo.rfo_root+1);
        rfo_root = optinfo.rfo_root;
      }
      close_PSIF();
    }

  /*for (i=0; i<dim+1; ++i) {
    if ((lambda[i] < 0.0) || (i <= rfo_root))
      if (fabs(rfo_mat[i][dim]) > 1.0e-4)
        for (j=0;j<dim+1;++j)
          rfo_mat[i][j] /= rfo_mat[i][dim];
  } */

    // print out lowest energy evects
    for (i=0; i<dim+1; ++i) {
      if ((lambda[i] < 0.0) || (i <rfo_root)) {
        fprintf(outfile,"RFO eigenvalue %d: %15.10lf (or 2*%-15.10lf)\n", i+1, lambda[i],lambda[i]/2);
        fprintf(outfile,"eigenvector:\n");
        print_mat2(&(rfo_mat[i]),1,dim+1,outfile);
      }
    }

/* alternative algorithm (H-lambdaI)x + g = 0
    dq = init_array(dim);
    double **H_test, **H_inv_test;
    H_test = block_matrix(dim,dim);
    for (i=0; i<dim; ++i)
      for (j=0; j<dim; ++j)
        H_test[i][j] = H[i][j];
    for (i=0; i<dim; ++i)
      H_test[i][i] = H[i][i] - lambda[0];
    H_inv_test = symm_matrix_invert(H_test,dim,0,0);
    mmult(H_inv_test,0,&f_q,1,&dq,1,dim,dim,1,0);
    fprintf(outfile,"dq solved by (H- lambda I) x + g = 0\n");
    print_mat2(&dq,1,dim,outfile);
    dot_arr(f_q, dq, dim, &tval);
    fprintf(outfile,"g^T x = lambda ? = %15.10lf\n", - tval);
    free(dq);
    free_block(H_test);
    free_block(H_inv_test);
*/

    // get norm |x| and unit vector in the step direction
    rfo_u = init_array(dim);
    for (j=0; j<dim; ++j)
      rfo_u[j] = rfo_mat[rfo_root][j];
    dot_arr(rfo_u, rfo_u, dim, &tval);
    rfo_xnorm = sqrt(tval);
    normalize(&rfo_u, 1, dim);

    // get gradient and hessian in step direction
    dot_arr(f_q, rfo_u, dim, &rfo_g);
    rfo_g *= -1; // gradient not force
    rfo_h = 0;
    for (i=0; i<dim; ++i) {
      dot_arr(H[i], rfo_u, dim, &tval);
      rfo_h += rfo_u[i] * tval;
    }

    dq = init_array(dim);
    for (j=0; j<dim; ++j)
      dq[j] = rfo_mat[rfo_root][j];
    rfo_eval = lambda[rfo_root];
    free(lambda);
    free_block(rfo_mat);

    // for debugging
    DE_new = rfo_energy(rfo_xnorm, rfo_g, rfo_h);
    fprintf(outfile,"DE_RFO(rfo_xnorm = %8.3e) = %20.10lf\n", rfo_xnorm, DE_new);
    //for (rfo_t=0; rfo_t<2*rfo_xnorm; rfo_t += (rfo_xnorm/10.0))
    //  fprintf(outfile,"DE_RFO(rfo_t = %8.3e) = %20.10lf\n", rfo_t, rfo_energy(rfo_t, rfo_g, rfo_h));

    fprintf(outfile,"Projected energy change by RFO approximation: %20.10lf\n", DE_new);
    open_PSIF();
    sprintf(value_string,"DE prediction %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &DE_new, sizeof(double));
    sprintf(value_string,"Unit step U %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) rfo_u, dim*sizeof(double));
    sprintf(value_string,"Xnorm %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &rfo_xnorm, sizeof(double)); 
    sprintf(value_string,"Scalar gradient %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &rfo_g, sizeof(double)); 
    sprintf(value_string,"Scalar hesian %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &rfo_h, sizeof(double));
    close_PSIF();

    free(rfo_u);
  } // end take RFO step
  } // take forward geometry step


  fprintf(outfile, "\n\tInternal Coordinate Update in Ang or Rad, aJ/Ang or aJ/Rad\n");
  fprintf(outfile, "\t----------------------------------------------------------\n");
  fprintf(outfile, "\t       Value         Force        Displacement  New Value\n");
  fprintf(outfile, "\t----------------------------------------------------------\n");
  for (i=0;i<dim;++i)
    fprintf(outfile,"\t%2d %13.8lf %13.8lf %13.8lf %13.8lf\n", i+1, q[i], f_q[i], dq[i], q[i]+dq[i]);
  fprintf(outfile, "\t----------------------------------------------------------\n");
  fprintf(outfile,"\n\t   MAX force: %15.10lf   RMS force: %15.10lf\n", max_force, rms_force);

  if (optinfo.analytic_interfragment) { // do interfragment steps analytically

    for (isalc=0; isalc<dim; ++isalc) {
     
      // find any interfragment coordinate sets
      simple_id = symm.get_simple(isalc,0);
      simples.locate_id(simple_id,&intco_type,&sub_index,&sub_index2);
      // only execute for single (non-salc) interfragment coordinates
      if ( (symm.get_length(isalc) != 1) || (intco_type != FRAG_TYPE) ) continue;
      // only execute this code once for each fragment pair
      // assume this one time is for the distance coordinate entry
      if (sub_index2 != 0) continue;
      // don't do analytic steps unless all reference points (and all coordinates) are present
      // because we may have to fix their values in orient_fragment even if we are not optmizing
      // with them
      if ( (simples.frag.get_A_P(sub_index) != 3) || (simples.frag.get_B_P(sub_index) != 3) )
        continue;

      for (cnt=0, I=0; I<6; ++I) {
        inter_q[I] = 0;
        if ( simples.frag.get_coord_on(sub_index,I) ) {
          inter_q[I] = dq[isalc + cnt];
          dq[isalc + cnt] = 0; // so that back-transformation code if used, doesn't change q again
          ++cnt;
        }
      }

      dot_arr(inter_q,inter_q,6,&tval); // check if there are any non-zero displacements
      if (fabs(tval) < 1.0e-15) continue;

      // calculate new coordinates
      for (cnt=0, I=0; I<6; ++I) {
        if ( simples.frag.get_coord_on(sub_index,I) )
          inter_q[I] += q[isalc + cnt++];
      }
      if (optinfo.frag_dist_rho)
        inter_q[0]  = 1.0 / inter_q[0];
      inter_q[0]  /= _bohr2angstroms;
      for (I=1; I<6; ++I) inter_q[I] *= 180.0/_pi; 

      // get geometries of fragments
      geom = carts.get_coord_2d();

      a = simples.frag.get_A_natom(sub_index);
      geom_A = block_matrix(a,3);
      for (atom=0;atom<a;++atom)
        for (xyz=0;xyz<3;++xyz)
          geom_A[atom][xyz] = geom[simples.frag.get_A_atom(sub_index,atom)][xyz];

      b = simples.frag.get_B_natom(sub_index);
      geom_B = block_matrix(b,3);
      for (atom=0;atom<b;++atom)
        for (xyz=0;xyz<3;++xyz)
          geom_B[atom][xyz] = geom[simples.frag.get_B_atom(sub_index,atom)][xyz];

      // get reference point information for fragments
      weight_A = block_matrix(3,a);
      for (cnt=0; cnt<simples.frag.get_A_P(sub_index); ++cnt)
        for (atom=0;atom<a;++atom)
          weight_A[cnt][atom] = simples.frag.get_A_weight(sub_index,cnt,atom);

      weight_B = block_matrix(3,b);
      for (cnt=0; cnt<simples.frag.get_B_P(sub_index); ++cnt)
        for (atom=0;atom<b;++atom)
          weight_B[cnt][atom] = simples.frag.get_B_weight(sub_index,cnt,atom);

      fprintf(outfile,"\nAnalytically doing interfragment displacements\n");
      // move fragment B and put result back into main geometry matrix
      orient_fragment(a, b, simples.frag.get_A_P(sub_index), simples.frag.get_B_P(sub_index),
        geom_A, geom_B, weight_A, weight_B, inter_q[0], inter_q[1], inter_q[2], inter_q[3],
        inter_q[4], inter_q[5], outfile);

      for (atom=0; atom<b; ++atom)
        for (xyz=0; xyz<3; ++xyz)
          geom[simples.frag.get_B_atom(sub_index,atom)][xyz] = geom_B[atom][xyz];

      symmetrize_geom(geom[0]);
      carts.set_coord_2d(geom);
      free_block(geom);
      free_block(geom_A); free_block(geom_B);
      free_block(weight_A); free_block(weight_B);
    }
    //fprintf(outfile,"\nNew Cartesian Geometry in a.u. after analytic interfragment steps\n");
    //carts.print(1, outfile,0, disp_label, 0);
  } // if analytic_interfragment

  // Do displacements with iterative backtransformation
  for (i=0;i<dim;++i)
    q[i] += dq[i];

  success = false;
  // check to see if all displacements are zero
  dot_arr(dq,dq,dim,&tval);
  if (tval < 1.0e-15) {
    success = true;
    carts.print(PSIF_CHKPT,outfile,0,NULL,0);
  }

  if (!success) {
    /* new_geom will overwrite dq so send a copy */
    dq_to_new_geom = init_array(dim);
    for (i=0;i<dim;++i)
      dq_to_new_geom[i] = dq[i];
  
    strcpy(disp_label,"New Cartesian Geometry in a.u.");
    success = new_geom(carts,simples,symm,dq_to_new_geom,32,0,disp_label,0,0,djunk);
  }
  
  int retry = 1;
  while (!success) {
    fprintf(outfile,"\tWarning: halving displacement size and reducing back-transformation \
 convergence criteria.\n");
    optinfo.bt_dx_conv *= 5.0;
    optinfo.bt_dq_conv *= 5.0;
    for (i=0;i<dim;++i) {
      dq[i] = dq[i] / 2.0;
      dq_to_new_geom[i] = dq[i];
    }
    success = new_geom(carts,simples,symm,dq_to_new_geom,32,0,disp_label,0,0,djunk);
    ++retry;
    if ((retry == 5) || (retry == 2 && do_line_search)) {
      fprintf(outfile,"Giving up - unable to back-transform to new cartesian coordinates.\n");
      fclose(outfile);
      exit(PSI_RETURN_FAILURE);
    }
  }

  // Modify predicted energe change if displacement size was reduced
  if (retry > 2) {
    dot_arr(f_q, dq, dim, &tval);
    DE_old = -tval; // g^T x
    for (i=0; i<dim; ++i) {
      dot_arr(H[i], dq, dim, &tval);
      DE_old += 0.5 * dq[i] * tval; // 1/2 x^T H x
    }
    DE_old /= _hartree2J * 1.0e18;
    fprintf(outfile,"\tNew projected energy change (minus any analytically taken steps\n");
    fprintf(outfile,"\talong interfragment coordinates: %15.10lf\n", DE_old);
    open_PSIF();
    sprintf(value_string,"DE prediction %d", optinfo.iteration);
    psio_write_entry(PSIF_OPTKING, value_string, (char *) &(DE_old), sizeof(double));
    close_PSIF();
  }

  free(q); free(dq); free(f_q);
  if (!do_line_search) {
    free(H);
    optinfo.iteration += 1;
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
    close_PSIF();
  }
  delete [] djunk;
  delete [] disp_label;
  return 0;
}


bool line_search(cartesians &carts, int dim, double *dq) {
  int i, dim_carts;
  double *u, xnorm, g, h, M, a, b, c, tval, tval2, t, *old_x;
  char value_string[20];
  double E, E_old, DE, DE_old;

  fprintf(outfile,"\tTrying smaller step in same direction from previous geometry.\n");
      
  E = carts.get_energy();
  dim_carts = 3*carts.get_natom();
      
  open_PSIF();
  sprintf(value_string,"Energy %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &E_old, sizeof(double));
  DE = E - E_old;
  fprintf(outfile,"\t\tActual DE   : %20.10lf\n", DE);

  sprintf(value_string,"DE prediction %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &DE_old, sizeof(double));
  fprintf(outfile,"\t\tPredicted DE: %20.10lf\n", DE_old);

  u = init_array(dim);
  sprintf(value_string,"Unit step U %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) u, dim*sizeof(double));
  //fprintf(outfile,"Unit step U %d\n", optinfo.iteration-1);
  //for (i=0; i<dim; ++i)
  //  fprintf(outfile,"%15.10lf\n", u[i]);

  sprintf(value_string,"Xnorm %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &xnorm, sizeof(double));
  fprintf(outfile,"\t\tXnorm %d     : %20.10lf\n", optinfo.iteration-1, xnorm);

  sprintf(value_string,"Scalar gradient %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &g, sizeof(double));
  fprintf(outfile,"\t\tScalar gradient %d: %15.10lf\n", optinfo.iteration-1, g);

  sprintf(value_string,"Scalar hesian %d", optinfo.iteration-1);
  psio_read_entry(PSIF_OPTKING, value_string, (char *) &h, sizeof(double));
  fprintf(outfile,"\t\tScalar hessian %d : %15.10lf\n", optinfo.iteration-1, h);

  if (optinfo.rfo)
    fprintf(outfile,"\t\tCheck energy: %20.15lf\n", rfo_energy(xnorm, g, h));
  else
    fprintf(outfile,"\t\tCheck energy: %20.15lf\n", nr_energy(xnorm, g, h));
      
  // sign of third derivative is + if energy obtained was too HIGH, and - if energy
  // if energy was too low, which we will only worry about for TS searches
  M = 6 * (DE-DE_old) / (xnorm*xnorm *xnorm) * (_hartree2J*1.0e18); //aJ/Ang^3
  fprintf(outfile,"\t\tApproximate 3rd derivative, M : %15.10lf\n",M);
      
  bool success = false;
  // reduce step size in 5% increments and test magnitude of cubic term
  for (t = xnorm; t > 0.01*xnorm; t -= 0.05*xnorm) {
    //fprintf(outfile,"delta abs(g+0.5ht) = %10.5e\n", optinfo.step_energy_limit_back*fabs(g+0.5*h*t));
    //fprintf(outfile,"1/6 Mt^2 = %10.5e\n", M*t*t/6);
    if (M*t*t/6 < optinfo.step_energy_limit_back*fabs(g+0.5*h*t)) {
      success = true;
      break;
    }
  }
  fprintf(outfile,"\t\tNew scalar with which to scale unit step, t = %15.10lf\n", t);
  fprintf(outfile,"\nSetting step to previous step reduced to %.0lf%% of original.\n\n", t/xnorm*100);

  if (success) {
    // put last geometry into cartesian object
    old_x = init_array(dim_carts);
    sprintf(value_string,"Cartesian Values %d", optinfo.iteration-1);
    psio_read_entry(PSIF_OPTKING, value_string, (char *) old_x, dim_carts*sizeof(double));
    carts.set_coord(old_x);
    fprintf(outfile,"Setting cartesian coordinates to prevous step.\n");
    carts.print(5, outfile, 0, NULL, 0);
    free(old_x);

    // reduced step is dq = u * rho_t;
    for (i=0; i<dim; ++i)
      dq[i] = t * u[i]; 

    //fprintf(outfile,"dq's\n");
    //for (i=0; i<dim; ++i)
    //  fprintf(outfile,"%15.10lf\n", dq[i]);
  }

  // what needs to be changed from 1st try from previous geometries?
  // energy, gradient, hessian, unit step, unchanged
  if (optinfo.rfo)
    DE_old = rfo_energy(t, g, h);
  else
    DE_old = nr_energy(t, g, h);

  sprintf(value_string,"DE prediction %d", optinfo.iteration-1);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &DE_old, sizeof(double));
  fprintf(outfile,"\tNew Projected DE %d: %20.10lf\n", optinfo.iteration-1,DE_old); 
   
  sprintf(value_string,"Xnorm %d", optinfo.iteration-1);
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &t, sizeof(double));
  
  close_PSIF();
  fflush(outfile);
  return success;
}

}} /* namespace psi::optking */

