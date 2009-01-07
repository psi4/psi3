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
extern void empirical_H(internals &simples, salc_set &symm, cartesians &carts);
void fconst_init(cartesians &carts, internals &simples, salc_set &symm);
extern void compute_zmat(cartesians &carts, int *unique_zvars);
extern void print_zmat(FILE *outfile, int *unique_zvars);
extern double **compute_H(internals &simples, salc_set &symm, double **P, cartesians &carts);

int opt_step(cartesians &carts, internals &simples, salc_set &symm) {

  int xyz, i,j,k,ii,a,b, dim, dim_carts, success,nbfgs, nsimples, constraint,cnt;
  double **B, **G, **G2, **G_inv, **H_inv, **temp_mat, **u, **P;
  double **C, **T, **T2, **T3, R_AB, theta_A, theta_B, tau, phi_A, phi_B;
  double *temp_arr2, *temp_arr, *masses, **geom, *forces;
  double *f, *f_q, *dq, *dq_to_new_geom, *q, tval, tval2, scale, temp, *djunk;
  char *disp_label, *wfn, value_string[30], force_string[30];

  // for analytic interfragment coordinates
  int isalc, I, simple_id, intco_type, sub_index, sub_index2, atom;
  double **geom_A, **geom_B, **weight_A, **weight_B;
  double inter_q[6];

  dim_carts = carts.get_natom()*3;
  djunk = new double[dim_carts];
  disp_label = new char[MAX_LINELENGTH];

  if (symm.get_num() == 0)
    punt("No symmetric internal coordinates to optimize.\n");

  ip_string("WFN", &(wfn),0);
  fprintf(outfile,"\nCurrent %s energy before step %20.10lf\n", wfn, carts.get_energy());
  free(wfn);
  fprintf(outfile,"\nTaking geometry step number %d\n",optinfo.iteration+1);

  //*** Build transformation matrices
  dq = init_array(symm.get_num());
  q = compute_q(simples,symm);

  // build G = BuB^t
  B = compute_B(simples,symm);
  G = compute_G(B,symm.get_num(),carts);

  //fprintf(outfile,"B matrix\n");
  //print_mat(B, symm.get_num(), dim_carts,outfile);

  // compute G_inv
  fprintf(outfile,"\nBuB^t ");
  G_inv = symm_matrix_invert(G,symm.get_num(),1,optinfo.redundant);

  // setup the masses matrix, u
  masses = carts.get_mass();
  u = mass_mat(masses);
  free(masses);
  //u = unit_mat(dim_carts);
  //fprintf(outfile,"Mass Matrix\n");
  //print_mat(u, dim_carts,dim_carts,outfile);

  // get forces array in cartesian coordinates, f, (in aJ/Ang)
  f = carts.get_forces();
  //fprintf(outfile,"cartesian forces in aJ/Ang\n");
  //print_mat2(&f, 1, dim_carts, outfile);

  // compute forces in internal coordinates, f_q = G_inv B u f
  f_q = init_array(symm.get_num());
  temp_arr = init_array(symm.get_num());
  temp_arr2 = init_array(dim_carts);

  mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
  mmult(B,0,&temp_arr2,1,&temp_arr,1,symm.get_num(),dim_carts,1,0);
  mmult(G_inv,0,&temp_arr,1,&f_q,1,symm.get_num(),symm.get_num(),1,0);
  free_block(B);

  // fprintf(outfile,"internal forces (G_inv B u f) in aJ/A\n");
  // print_mat2(&f_q, 1, symm.get_num(), outfile);
  // for (i=0;i<symm.get_num();++i)
  //   f_q[i] = f_q[i] * sqrt(2) * _hartree2J * 1.0E18 / _bohr2angstroms;
  // fprintf(outfile, "%d %15.10lf\n", i, f_q[i] * _hartree2J / _bohr2angstroms * 1.0E18) ;

  // test by transforming f_q back to cartesian forces and compare
  // fprintf(outfile,"computed forces in cartesian coordinates aJ/Ang\n");
  // mmult(B,1,&f_q,1,&temp_arr2,1,dim_carts,symm.get_num(),1,0);
  // print_mat2(&temp_arr2, 1, dim_carts, outfile);

  free(f);
  free(temp_arr);
  free(temp_arr2);
  free_block(u);

  // Setup projection matrix P = G * G-
  // for inversion of Hessian with redundant internals
  P = block_matrix(symm.get_num(),symm.get_num());
  mmult(G,0,G_inv,0,P,0,symm.get_num(),symm.get_num(),symm.get_num(),0); 
  free_block(G);
  free_block(G_inv);
  
  // Now apply constraints
  if (optinfo.constraints_present) {
    dim = symm.get_num();
    // C has 1's on diagonal for constraints and 0's elsewhere
    C = block_matrix(dim,dim);
    if (optinfo.redundant) {
      for (i=0;i<optinfo.nconstraints;++i)
        C[optinfo.constraints[i]][optinfo.constraints[i]] = 1.0;
    }
    else if (optinfo.delocalize) {

      double i_overlap, j_overlap;
      for (b=0; b<optinfo.nconstraints; ++b) {
        constraint = simples.index_to_id(optinfo.constraints[b]);

        for (i=0; i<dim; ++i)
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
     // print_mat2(C,dim,dim,outfile);
     // fflush(outfile);

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

   // Project forces f_q = P f_q'
   temp_arr = init_array(symm.get_num());
   mmult(P,0,&f_q,1,&temp_arr,1,symm.get_num(), symm.get_num(),1,0);
   for(i=0; i<symm.get_num(); ++i)
     f_q[i] = temp_arr[i];
   free(temp_arr);

  // Write Values and Forces of internals to opt.aux for later
  open_PSIF();
  nbfgs = 0;
  if (psio_tocscan(PSIF_OPTKING, "Num. of Previous Entries") != NULL)
    psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &nbfgs, sizeof(int));

  sprintf(value_string,"Previous Internal Values %d", nbfgs);
  sprintf(force_string,"Previous Internal Forces %d", nbfgs);
       
  psio_write_entry(PSIF_OPTKING, value_string, (char *) &(q[0]),
      symm.get_num()* sizeof(double));
  psio_write_entry(PSIF_OPTKING, force_string, (char *) &(f_q[0]),
      symm.get_num()* sizeof(double));

  double *coord;
  coord = carts.get_coord();

  sprintf(value_string,"Previous Cartesian Values %d", nbfgs);
  psio_write_entry(PSIF_OPTKING, value_string,
      (char *) &(coord[0]), 3*carts.get_natom()*sizeof(double));
  free(coord);

  ++nbfgs;
  psio_write_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &nbfgs, sizeof(int));
  close_PSIF();

  // make sure some force constant are in PSIF
  fconst_init(carts, simples, symm);

  // Read in Hessian and update it if necessary from opt.aux
  H_inv = compute_H(simples,symm,P,carts);
  free_block(P);

  // Computing internal coordinate displacements H_inv f = dq, and new q
  mmult(H_inv,0,&f_q,1,&dq,1,symm.get_num(),symm.get_num(),1,0);
  free_block(H_inv);

  // Fix any interfragment torsions that pass through 180.0
/*
  int simple, intco_type, sub_index, J, simple_b;
  int A_natom, B_natom, *A_atom, *B_atom;
  for (i=0;i<symm.get_num();++i) {
    if (symm.get_length(i) == 1) { // it's a simple
      simple = symm.get_simple(i,0);
      simples.locate_id(simple,&intco_type,&sub_index);
      if (intco_type == FRAG_TYPE) { // it's a fragment coordinate
        J = simples.frag.get_J(sub_index);
        A_natom = simples.frag.get_A_natom(sub_index);
        B_natom = simples.frag.get_B_natom(sub_index);
        A_atom = init_int_array(A_natom);
        B_atom = init_int_array(B_natom);
        for (a=0; a<A_natom; ++a)
          A_atom[a] = simples.frag.get_A_atom(sub_index,a);
        for (b=0; b<B_natom; ++b)
          B_atom[b] = simples.frag.get_B_atom(sub_index,b);
        if (J==2) { // polar angle on B is out-of-range 0-pi
          if ( (q[i]+dq[i] > _pi) || (q[i]+dq[i]<0.0)) {
            fprintf(outfile,"Polar angle on B is out-of-range, switching related torsions.\n");
            dq[i] = 0.0; // leave bond angle alone, switch the connected dihedrals
            simple_b = simples.frag.get_id_from_atoms(A_natom,B_natom,A_atom,B_atom,3);
            for (ii=0;ii<symm.get_num();++ii) {
              if (simple_b == symm.get_simple(ii,0)) {
                q[ii] *= -1.0;
                dq[ii] = 0.0;
              }
            }
            simple_b = simples.frag.get_id_from_atoms(A_natom,B_natom,A_atom,B_atom,5);
            for (ii=0;ii<symm.get_num();++ii) {
              if (simple_b == symm.get_simple(ii,0)) {
                q[ii] *= -1.0;
                dq[ii] = 0.0;
              }
            }
          }
        } // polar angle B
        free(A_atom); free(B_atom);
      } //fragment type
    }
  }
*/
  // Test forces to see if geometry is optimized
  tval = 0.0; 
  tval2 = fabs(f_q[0]);
  for (i=0;i<symm.get_num();++i) {
    tval += SQR(f_q[i]);
    if (fabs(f_q[i]) > tval2) tval2 = fabs(f_q[i]);
  } 
  tval = tval/((double) symm.get_num());
  tval = sqrt(tval);
      
  if (tval2 < optinfo.conv) {
    fprintf(outfile, "\nInternal Coordinates and Forces in Ang or Rad, aJ/Ang or aJ/Rad\n");
    fprintf(outfile, "       Value         Force        \n");
    for (i=0;i<symm.get_num();++i)
      fprintf(outfile,"%2d %13.8lf %13.8lf \n", i+1,  q[i],  f_q[i] );
    fprintf(outfile,"\nMAX force is < %5.1e.  Optimization is complete.\n",
            optinfo.conv);
    ip_string("WFN", &(wfn),0);
    fprintf(outfile,"Final %s energy is %15.10lf\n", wfn, carts.get_energy());
    fprintf(stderr,"\n  OPTKING:  optimization is complete\n");
    fprintf(outfile,"The Optimized geometry in a.u.\n");
    carts.print(12,outfile,0,disp_label,0);
    fprintf(outfile,"\nThe Optimized geometry in Angstrom\n");
    carts.print(13,outfile,0,disp_label,0);
    fprintf(outfile, "\n");
    fflush(outfile);

    if (optinfo.zmat) {
      int *unique_zvars;
      unique_zvars = init_int_array(MAX_ZVARS);
      compute_zmat(carts, unique_zvars);
      print_zmat(outfile, unique_zvars);
      free(unique_zvars);
      fprintf(outfile,"\n");
    }
    free(wfn); free(f_q); free(q); free(dq);
    fprintf(stderr,"\n  Returning code %d\n", PSI_RETURN_ENDLOOP);
    return(PSI_RETURN_ENDLOOP);
  } /* end converged geometry */

  // for now, lets not scale steps for interfragment only optimizations
  if (!optinfo.freeze_intrafragment) {
    /* determine scale factor needed to keep step less than 10% of q if q big
       or less than 0.1 if q not so big, hopefully better solution coming soon
       double cut = STEP_LIMIT / STEP_PERCENT; */
    scale = 1.0;
    for (i=0;i<symm.get_num();++i) {
      if (fabs(q[i]) > STEP_LIMIT) { /* intco is not very small */
        if (fabs(dq[i]) > STEP_PERCENT*fabs(q[i])) { /* step is too large */
          temp = STEP_PERCENT*fabs(q[i])/fabs(dq[i]);
        }
        else
          temp = 1;
      }
      else { /* intco is very small */
        if (fabs(dq[i]) < STEP_LIMIT)
          temp = 1.0; /* step is small enough */
        else
          temp = STEP_LIMIT / fabs(dq[i]);
      }
      if (temp < scale){
        scale = temp;
      }
    } 
    fprintf(outfile,"\nScaling displacements by %lf\n",scale);
    for (i=0;i<symm.get_num();++i) {
      dq[i] = dq[i] * scale;   
    }   
  }

  // avoid piddly symmetry breaking
  for (i=0;i<symm.get_num();++i)
    if (fabs(dq[i]) < MIN_DQ_STEP) dq[i] = 0.0;

  fprintf(outfile, "\nInternal Coordinate Update in Ang or Rad, aJ/Ang or aJ/Rad\n");
  fprintf(outfile, "       Value         Force        Displacement  New Value\n");
  for (i=0;i<symm.get_num();++i)
    fprintf(outfile,"%2d %13.8lf %13.8lf %13.8lf %13.8lf\n",
        i+1,  q[i],  f_q[i],  dq[i],  q[i]+dq[i]);
  fprintf(outfile,"   MAX force: %15.10lf   RMS force: %15.10lf\n",tval2,tval);
  free(f_q);

  printf("Energy: %15.10lf MAX force: %6.2e RMS force: %6.2e\n",carts.get_energy(),tval2,tval);

  if (optinfo.analytic_interfragment) { // do interfragment steps analytically

    for (isalc=0; isalc<symm.get_num(); ++isalc) {
     
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

      // limit step sizes some more ?
/*
      if (!optinfo.frag_dist_rho && (fabs(inter_q[0] > 0.2))) {
        fprintf(outfile,"Limiting change in coordinate %d to 0.2 angstroms.\n",simple_id);
        inter_q[0] = 0.2 * inter_q[0] / fabs(inter_q[0]);
      }
      for (I=1;I<6;++I) {
        if ( fabs(inter_q[I])*180.0/_pi > 10.0) {
          fprintf(outfile,"Limiting change in coordinate %d to 10 degrees.\n",simple_id);
          inter_q[I] = 10.0/180.0*_pi * inter_q[I] / fabs(inter_q[I]);
        }
      }
*/

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
  if (!optinfo.freeze_intrafragment) {

    for (i=0;i<symm.get_num();++i)
      q[i] += dq[i];

/* fprintf(outfile,"dq values for iterative back-trans\n");
for (i=0;i<symm.get_num();++i) fprintf(outfile,"%15.10lf\n",dq[i]); */

    /* new_geom will overwrite dq so send a copy */
    dq_to_new_geom = init_array(symm.get_num());
    for (i=0;i<symm.get_num();++i)
      dq_to_new_geom[i] = dq[i];
  
    strcpy(disp_label,"New Cartesian Geometry in a.u.");
    success = new_geom(carts,simples,symm,dq_to_new_geom,32,0,disp_label,0,0,djunk);
  
    int retry = 1;
    while (!success) {
      fprintf(outfile,"\tWarning: halving displacement size and reducing back-transformation \
 convergence criteria.\n");
      optinfo.bt_dx_conv *= 5.0;
      optinfo.bt_dq_conv *= 5.0;
      for (i=0;i<symm.get_num();++i) {
        dq[i] = dq[i] / 2.0;
        dq_to_new_geom[i] = dq[i];
      }
      success = new_geom(carts,simples,symm,dq_to_new_geom,32,0,disp_label,0,0,djunk);
      ++retry;
      if (retry == 15) {
        fprintf(outfile,"Giving up - unable to back-transform to new cartesian coordinates.\n");
        fclose(outfile);
        exit(PSI_RETURN_FAILURE);
      }
    }
  } // use iterative back-transformation

  free(q); free(dq);
  optinfo.iteration += 1;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Iteration", (char *) &(optinfo.iteration),sizeof(int));
  close_PSIF();
  delete [] djunk;
  delete [] disp_label;
  return 0;
}

}} /* namespace psi::optking */

