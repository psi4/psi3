/*! \file 
    \ingroup (ORIENT_FRAGMENTS)
    \brief arrange multiple fragments based on interfragment coordinates 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <cmath>
#include "input.h"
#include "global.h"
#include "defines.h"
#include "physconst.h"

namespace psi { namespace input {

void next_point(double *A, double *B, double *C, double R_CD, double theta_BCD, double phi_ABCD, double *D);
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);

void orient_fragments()
{
  int i, j, errcod, f, ival, pts, xyz, frag_dist_rho=0;
  int *P; // number of reference points in each fragment
  char *ref_pts_lbl, error_message[80], frag_lbl[80];
  double ***ref_pts_lc, B_angle, R_B1B2, R_B2B3, e12[3], e12b[3], e12c[3], e12d[3], erot[3];
  double R_AB, theta_A, theta_B, phi_A, phi_B, tau, tval, norm, **ref_B_final, **B;
  double **ref_A, **ref_B, **c, sign, cross1[3], cross2[3], cross3[3], phi2, phi3;

  ip_boolean("FRAGMENT_DISTANCE_INVERSE", &(frag_dist_rho),0);

  /* a reference point for a fragment is a linear combination of atom positions - each
    fragment needs 3 reference points, unless the fragment has 1 atoms, 2 atoms, or is linear */
  ref_pts_lc = (double ***) malloc(nfragments*sizeof(double **));
  P = (int *) malloc(nfragments*sizeof(int));
  for (f=0; f<nfragments; ++f) {
    if (frag_num_atoms[f]==1)      P[f] = 1;
    else if (frag_num_atoms[f]==2) P[f] = 2;
    else P[f] = 3;
    ref_pts_lc[f] = block_matrix(P[f],frag_num_atoms[f]);
  }

  ref_pts_lbl = (char *) malloc(30*sizeof(char)); 
  for (f=0; f<nfragments; ++f) {
    if (f == 0)
      sprintf(ref_pts_lbl,"GEOMETRY_REF_POINTS");
    else
      sprintf(ref_pts_lbl,"GEOMETRY%d_REF_POINTS",f+1);

    /* Read in reference points linear combinations - normalize them if necessary */
    errcod = ip_count(ref_pts_lbl, &ival, 0);
    if (ival < P[f]) {
      fprintf(outfile,"Caution: Too few reference pts for full geometry specification, unless fragment %d is non-linear\n",f+1);
      P[f] = ival;
    }
    else if (ival > P[f]) {
      sprintf(error_message,"%s requires %d reference points!", ref_pts_lbl, P[f]);
      punt(error_message);
    }
    for (i=0; i<P[f]; ++i) {
      errcod = ip_count(ref_pts_lbl, &ival, 1, i);
      if (ival != frag_num_atoms[f]) {
        sprintf(error_message,"Fragment: %d, Reference point: %d requires %d entries!", f+1, i+1, frag_num_atoms[f]);
        punt(error_message);
      }
      for (j=0; j<frag_num_atoms[f]; ++j)
        ip_data(ref_pts_lbl, "%lf", &ref_pts_lc[f][i][j], 2, i, j);
    }
    for (i=0; i<P[f]; ++i) {
      tval = 0.0;
      for (j=0; j<frag_num_atoms[f]; ++j)
        tval += fabs(ref_pts_lc[f][i][j]);
      tval = 1.0/tval;
      for (j=0; j<frag_num_atoms[f]; ++j)
        ref_pts_lc[f][i][j] *= tval;
    }
    fprintf(outfile,"Linear combinations which specify reference points for fragment %d\n", f);
    print_mat(ref_pts_lc[f],P[f],frag_num_atoms[f],outfile);
  }
 
  ref_A = block_matrix(3,3);
  ref_B = block_matrix(3,3);
  ref_B_final = block_matrix(3,3);

  for (f=1; f<nfragments; ++f) { /* leave first fragment alone */

    for (i=0; i<3; ++i)
      for (xyz=0; xyz<3; ++xyz) {
        ref_A[i][xyz] = 0.0;
        ref_B[i][xyz] = 0.0;
      }

    /* compute xyz coordinates for reference points on fragment A and B */
    for (pts=0; pts<P[f-1]; ++pts)
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<frag_num_atoms[f-1];++i)
          ref_A[pts][xyz] += ref_pts_lc[f-1][pts][i] * geometry[frag_atom[f-1]+i][xyz];
    /* add missing ref pts anywhere - just as long as they are not collinear, we're OK */
    if (P[f-1] < 3) {
      for (xyz=0; xyz<3; ++xyz)
        ref_A[2][xyz] = (xyz+1)/_pi;
    }
    if (P[f-1] < 2) {
      for (xyz=0; xyz<3; ++xyz)
        ref_A[1][xyz] = (xyz+1)/(2*_pi);
    }

    B = block_matrix(frag_num_atoms[f],3);
    for (i=0; i<frag_num_atoms[f];++i)
      for (xyz=0; xyz<3; ++xyz)
        B[i][xyz] = geometry[frag_atom[f]+i][xyz];

fprintf(outfile,"Coordinates for fragment:\n");
print_mat(B,frag_num_atoms[f],3,outfile);

    for (pts=0; pts<P[f]; ++pts)
      for (xyz=0; xyz<3; ++xyz)
        for (i=0; i<frag_num_atoms[f];++i)
          ref_B[pts][xyz] += ref_pts_lc[f][pts][i] * B[i][xyz];

    fprintf(outfile,"Coordinates for reference points on fragment A\n");
    print_mat(ref_A,P[f],3,outfile);
    fprintf(outfile,"Coordinates for reference points on fragment B (original) \n");
    print_mat(ref_B,P[f],3,outfile);

    /* read in interfragment coordinates A->B; B will be moved to fit these */
    sprintf(frag_lbl, "FRAGMENTS%d%d", f, f+1);
    if (!ip_exist(frag_lbl,0))
      punt("You need to specify interfragment coordinates.");
    ip_count(frag_lbl, &ival, 0);
    if (ival != 6)
      punt("Interfragment coordinates should be given in sets of 6 (add 0's if necessary)");
    ip_data(frag_lbl, "%lf", &R_AB   , 1, 0);
    if (frag_dist_rho)
      R_AB = 1.0/R_AB * conv_factor; /* convert to au - allow 1/R coordinate */
    else 
      R_AB = R_AB * conv_factor;
    ip_data(frag_lbl, "%lf", &theta_A, 1, 1);
    ip_data(frag_lbl, "%lf", &theta_B, 1, 2);
    ip_data(frag_lbl, "%lf", &tau    , 1, 3);
    ip_data(frag_lbl, "%lf", &phi_A  , 1, 4);
    ip_data(frag_lbl, "%lf", &phi_B  , 1, 5);
    fprintf(outfile,"\nGiven %d-%d interfragment coordinates:\n", f, f+1);
    fprintf(outfile,"\t(1/)R_AB:%10.5f, theta_A:%10.5f, theta_B:%10.5f\n", R_AB, theta_A, theta_B);
    fprintf(outfile,"\t     tau:%10.5f,   phi_A:%10.5f,   phi_B:%10.5f\n", tau, phi_A, phi_B);

    /* compute B1-B2 distance, B2-B3 distance, and B1-B2-B3 angle */
    R_B1B2 = 0.0;
    if (P[f]>1) {
      for (xyz=0; xyz<3; ++xyz)
        R_B1B2 += (ref_B[1][xyz]-ref_B[0][xyz])*(ref_B[1][xyz]-ref_B[0][xyz]);
      R_B1B2 = sqrt(R_B1B2);
    }
    R_B2B3 = 0.0;
    B_angle = 0.0;
    if (P[f]>2) {
      for (xyz=0; xyz<3; ++xyz)
        R_B2B3 += (ref_B[2][xyz]-ref_B[1][xyz])*(ref_B[2][xyz]-ref_B[1][xyz]);
      R_B2B3 = sqrt(R_B2B3);
      unit_vec(ref_B[1],ref_B[0],e12);
      unit_vec(ref_B[1],ref_B[2],e12b);
      B_angle = acos(dot_prod(e12,e12b))*180.0/_pi;
    }

    /* determine location of reference pts for B in coordinate system of A */
    next_point(ref_A[2], ref_A[1], ref_A[0], R_AB, theta_A, phi_A, ref_B_final[0]);
    if (P[f]>1)
      next_point(ref_A[1], ref_A[0], ref_B_final[0], R_B1B2, theta_B, tau, ref_B_final[1]);
    if (P[f]>2)
      next_point(ref_A[0], ref_B_final[0], ref_B_final[1], R_B2B3, B_angle, phi_B, ref_B_final[2]);

    fprintf(outfile,"Target reference points for fragment %d\n", f+1);
    print_mat(ref_B_final,P[f],3,outfile);

    /* translate ref_B to place B1 in correct location */
    for (xyz=0; xyz<3; ++xyz) {
      tval = ref_B_final[0][xyz] - ref_B[0][xyz];
      for (i=0; i<frag_num_atoms[f]; ++i)
        B[i][xyz] += tval;
    }

fprintf(outfile,"Coordinates for fragment:\n");
print_mat(B,frag_num_atoms[f],3,outfile);

    for (pts=0; pts<P[f]; ++pts)
      for (xyz=0; xyz<3; ++xyz) {
        ref_B[pts][xyz] = 0.0;
        for (i=0; i<frag_num_atoms[f];++i)
          ref_B[pts][xyz] += ref_pts_lc[f][pts][i] * B[i][xyz];
      }

    fprintf(outfile,"Reference points after translation (to fix point B1):\n");
    print_mat(ref_B,P[f],3,outfile);

    if (P[f]>1) { /* move fragment B to place reference point B2 in correct location */
      /* Determine rotational angle and axis */
      unit_vec(ref_B[1],       ref_B[0], e12);  /* v B1->B2 */
      unit_vec(ref_B_final[1], ref_B[0], e12b); /* v B1->B2_final */
      B_angle = acos(dot_prod(e12b,e12));
      fprintf(outfile,"Rotation by %f degrees (to fix point B2)\n", 180.0*B_angle/_pi); 
      if (fabs(B_angle) > 1.0e-14) { 
        cross_prod(e12,e12b,erot);

        /* Move B to put B1 at origin */
        for (xyz=0; xyz<3; ++xyz) 
          for (i=0; i<frag_num_atoms[f];++i)
            B[i][xyz] -= ref_B[0][xyz];
  
        /* Rotate B */
        rotate_3D(erot, B_angle, B, frag_num_atoms[f]);
  
        /* Move B back to coordinate system of A */
        for (xyz=0; xyz<3; ++xyz) 
          for (i=0; i<frag_num_atoms[f];++i)
            B[i][xyz] += ref_B[0][xyz];
  
  fprintf(outfile,"Coordinates for fragment:\n");
  print_mat(B,frag_num_atoms[f],3,outfile);
      
        /* Check location of reference points now */
        for (pts=0; pts<P[f]; ++pts)
          for (xyz=0; xyz<3; ++xyz) {
            ref_B[pts][xyz] = 0.0;
            for (i=0; i<frag_num_atoms[f];++i)
              ref_B[pts][xyz] += ref_pts_lc[f][pts][i] * B[i][xyz];
          }

        fprintf(outfile,"Reference points after rotation (to fix point B2) \n");
        print_mat(ref_B,P[f],3,outfile);
      }
    }

    if (P[f]==3) { /* move fragment B to place reference point B3 in correct location */
      /* Determine rotational angle and axis */
      unit_vec(ref_B[1], ref_B[0], erot);  /* B1 -> B2 is rotation axis */

      /* Calculate B3-B1-B2-B3' torsion angle */
      unit_vec(ref_B[2], ref_B[0], e12);  /* v B1->B3 */
      unit_vec(ref_B[1], ref_B[0], e12b); /* v B1->B2 */
      phi2 = acos(dot_prod(e12,e12b));
      unit_vec(ref_B[0], ref_B[1], e12c);  /* v B2->B1 */
      unit_vec(ref_B_final[2], ref_B[1], e12d); /* v B2->B3' */
      phi3 = acos(dot_prod(e12c,e12d));

      cross_prod(e12 , e12b, cross1) ; /* B3->B1 x B1->B2 */
      cross_prod(e12c, e12d, cross2) ; /* B1->B2 x B2->B3 */
      tval = dot_prod(cross1, cross2) ;

      if ((sin(phi2) > 0.00001) && (sin(phi3) > 0.00001)) {
        tval /= sin(phi2) ;
        tval /= sin(phi3) ;
      }
      else tval = 2.0;

      if (tval > 0.99999) B_angle = 0.0000;
      else if (tval < -0.99999) B_angle = _pi;
      else B_angle = acos(tval) ;

      sign = 1.0; /* check sign */
      cross_prod(cross1, cross2, cross3);
      norm = sqrt(dot_prod(cross3, cross3));
      if (fabs(norm) > 0.00001) {
        for (xyz=0; xyz<3; ++xyz)
          cross3[xyz] *= 1.0/norm;
        tval = dot_prod(cross3, e12b);
        if (tval < 0.0) sign = -1.0;
      }
      B_angle *= sign;

      fprintf(outfile,"Rotation by %f degrees (to fix point B3)\n", 180.0*B_angle/_pi); 

      /* Move B to put B2 at origin */
      for (xyz=0; xyz<3; ++xyz) 
        for (i=0; i<frag_num_atoms[f];++i)
          B[i][xyz] -= ref_B[1][xyz];

      rotate_3D(erot, B_angle, B, frag_num_atoms[f]);

      /* Translate B1 back to coordinate system of A */
      for (xyz=0; xyz<3; ++xyz) 
        for (i=0; i<frag_num_atoms[f];++i)
          B[i][xyz] += ref_B[1][xyz];

fprintf(outfile,"Coordinates for fragment:\n");
print_mat(B,frag_num_atoms[f],3,outfile);
    
      for (pts=0; pts<P[f]; ++pts)
        for (xyz=0; xyz<3; ++xyz) {
          ref_B[pts][xyz] = 0.0;
          for (i=0; i<frag_num_atoms[f];++i)
            ref_B[pts][xyz] += ref_pts_lc[f][pts][i] * B[i][xyz];
        }

      fprintf(outfile,"Reference points on B after rotation for B2 \n");
      print_mat(ref_B,P[f],3,outfile);
    }

    /* check to see if desired reference points were obtained */
    tval = 0.0;
    for (i=0; i<P[f]; ++i)
      for (xyz=0; xyz<3; ++xyz)
        tval += fabs(ref_B[i][xyz] - ref_B_final[i][xyz]);
    if (tval > 1.0e10)
      punt("Unable to construct starting geometry.");
    else
      fprintf(outfile,"Starting interfragment coordinates achieved.\n");

    for (i=0;i<frag_num_atoms[f];++i)
      for (xyz=0; xyz<3; ++xyz)
        geometry[frag_atom[f]+i][xyz] = B[i][xyz];

    free_block(B);
  }

  free_block(ref_A);
  free_block(ref_B);
  free(ref_pts_lbl);
  free(P);
  for (f=0; f<nfragments; ++f)
    free_block(ref_pts_lc[f]);
  free(ref_pts_lc);
  return;
}

/* Given xyz coordinates for three points and R, theta, and phi, returns the
coordinates a fourth point; angles should enter function in degrees */

void next_point(double *A, double *B, double *C, double R_CD, double theta_BCD,
  double phi_ABCD, double *D)
{
  double eAB[3],eBC[3],eX[3],eY[3], cosABC, sinABC;
  int xyz;

  theta_BCD *= _pi/180.0;
  phi_ABCD *= _pi/180.0;

  unit_vec(B,A,eAB); /* vector B->A */
  unit_vec(C,B,eBC); /* vector C->B */
  cosABC = -dot_prod(eBC,eAB);

  sinABC = sqrt(1 - (cosABC * cosABC) );
  if ( (sinABC - LINEAR_CUTOFF) < 0.0 ) {
    fprintf(outfile,"  Reference points cannot be collinear.");
    punt("Invalid interfragment coordinates");
  }

  cross_prod(eAB,eBC,eY);
  for(xyz=0;xyz<3;xyz++)
    eY[xyz] /= sinABC;
  cross_prod(eY,eBC,eX);
  for (xyz=0;xyz<3;xyz++)
    D[xyz] = C[xyz] + R_CD * ( - eBC[xyz] * cos(theta_BCD) +
                                 eX[xyz] * sin(theta_BCD) * cos(phi_ABCD) +
                                 eY[xyz] * sin(theta_BCD) * sin(phi_ABCD) );
  return;
}

/*** SYM_MATRIX_INVERT inverts a matrix by diagonalization
 *  **A = matrix to be inverted
 *  dim = dimension of A
 *  print_det = print determinant if 1, nothing if 0
 *  redundant = zero eigenvalues allowed if 1
 *  returns: inverse of A ***/
#define REDUNDANT_EVAL_TOL (1.0e-10)
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

  A_inv   = block_matrix(dim,dim);
  A_temp  = block_matrix(dim,dim);
  A_vects = block_matrix(dim,dim);
  A_vals  = init_array(dim);

  sq_rsp(dim,dim,A,A_vals,1,A_vects,1.0e-14);

  if (redundant == 0) {
    for (i=0;i<dim;++i) { 
      det *= A_vals[i];
      A_inv[i][i] = 1.0/A_vals[i];
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
    if (fabs(det) < 1E-10)
      fprintf(outfile,"Warning, determinant is small: %10.6e\n",det);
  }
  else {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
      if (fabs(A_vals[i]) > REDUNDANT_EVAL_TOL)
        A_inv[i][i] = 1.0/A_vals[i];
      else
        A_inv[i][i] = 0.0;
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
  } 
      
  mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);
    
  free(A_vals);
  free_block(A_vects);
  free_block(A_temp);
  return A_inv;
}

}} // namespace psi::input
