/* six coordinates for non-bonded fragments */

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <physconst.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "fragment.h"

void fragment_class::print(FILE *fp_out, int print_values, int print_weights) {
  int i,a,b;
  fprintf(fp_out,"    (%d ", id);
  fprintf(fp_out," %d ", J+1);
  fprintf(fp_out,"( ", id);
  for (i=0; i<A_natom; ++i)
    fprintf(fp_out,"%d ", A_atom[i]+1);
  fprintf(outfile,")");
  fprintf(fp_out,"   ( ", id);
  for (i=0; i<B_natom; ++i)
    fprintf(fp_out,"%d ", B_atom[i]+1);
  fprintf(outfile,") )\n");

  if (print_weights) {
    fprintf(outfile,"    Weights for fragment A reference points **\n");
    for (i=0; i<3; ++i) {
        fprintf(outfile,"     (");
      for (a=0; a<A_natom; ++a)
        fprintf(outfile," %.8lf", A_weight[i*A_natom+a]);
      fprintf(outfile,")\n");
    }
    fprintf(outfile,"    Weights for fragment B reference points **\n");
    for (i=0; i<3; ++i) {
        fprintf(outfile,"     (");
      for (b=0; b<B_natom; ++b)
        fprintf(outfile," %.8lf", B_weight[i*B_natom+b]);
      fprintf(outfile,")\n");
    }
  }

  if (print_values) {
    if (J==0) {
      if (optinfo.frag_dist_rho)
        fprintf(fp_out, "   J=1, 1/R(AB)  = %12.6lf\n", value);
      else
        fprintf(fp_out, "   J=1, R(AB)    = %12.6lf\n", value);
    }
    else if (J==1)
      fprintf(fp_out, "   J=2, theta-A  = %12.6lf\n", value);
    else if (J==2)
      fprintf(fp_out, "   J=3, theta-B  = %12.6lf\n", value);
    else if (J==3)
      fprintf(fp_out, "   J=4, tau A-B  = %12.6lf\n", value);
    else if (J==4)
      fprintf(fp_out, "   J=5, chi-A    = %12.6lf\n", value);
    else if (J==5)
      fprintf(fp_out, "   J=6, chi-B    = %12.6lf\n", value);
  }
  else fprintf(outfile,"\n");
}

double fragment_class::get_val_A_or_rad(void) {
  double tval;
  if (J==0)
    tval = value;
  else if (J==1)
    tval = value *_pi/180.0;
  else if (J==2)
    tval = value *_pi/180.0;
  else if (J==3)
    tval = value *_pi/180.0;
  else if (J==4)
    tval = value *_pi/180.0;
  else if (J==5)
    tval = value *_pi/180.0;
  return tval;
}

void fragment_class::compute(double *geom) {
  int xyz, k, i;
  double **dkA, **dkB; /* location of reference points */
  double *e12A, *e12B, *e32A, *e32B, *eRA, *eRB, *v, *v2, *v3;
  double  d12A,  d12B,  d32A,  d32B, R;
  double dot, alpha_A, alpha_B, theta_A, theta_B;

  dkA = block_matrix(3,3);
  dkB = block_matrix(3,3);
  e12A = init_array(3);
  e12B = init_array(3);
  e32A = init_array(3);
  e32B = init_array(3);
  eRA = init_array(3);
  eRB = init_array(3);
  v = init_array(3);
  v2 = init_array(3);
  v3 = init_array(3);

  /* compute reference points within each fragment */
  for (k=0; k<3; ++k) { /* k = point 1, 2 or 3 */
    for (xyz=0; xyz<3; ++xyz) {
      for (i=0; i<A_natom; ++i) 
        dkA[k][xyz] += A_weight[k*A_natom+i] * geom[3*A_atom[i]+xyz];
      for (i=0; i<B_natom; ++i) 
        dkB[k][xyz] += B_weight[k*B_natom+i] * geom[3*B_atom[i]+xyz];
    }
  }
//print_mat(dkA, 3, 3, outfile);
//print_mat(dkB, 3, 3, outfile);

  for (xyz=0; xyz<3; ++xyz) {
    e12A[xyz] = dkA[1][xyz] - dkA[0][xyz];
    e12B[xyz] = dkB[1][xyz] - dkB[0][xyz];
    e32A[xyz] = dkA[1][xyz] - dkA[2][xyz];
    e32B[xyz] = dkB[1][xyz] - dkB[2][xyz];
    eRA[xyz]  = dkB[0][xyz] - dkA[0][xyz]; 
    eRB[xyz]  = dkA[0][xyz] - dkB[0][xyz];
  }
  d12A = sqrt( SQR(e12A[0]) +SQR(e12A[1]) +SQR(e12A[2]) );
  d12B = sqrt( SQR(e12B[0]) +SQR(e12B[1]) +SQR(e12B[2]) );
  d32A = sqrt( SQR(e32A[0]) +SQR(e32A[1]) +SQR(e32A[2]) );
  d32B = sqrt( SQR(e32B[0]) +SQR(e32B[1]) +SQR(e32B[2]) );
  R    = sqrt( SQR( eRA[0]) +SQR( eRA[1]) +SQR( eRA[2]) );

  scalar_mult(1/d12A, e12A, 3);
  scalar_mult(1/d12B, e12B, 3);
  scalar_mult(1/d32A, e32A, 3);
  scalar_mult(1/d32B, e32B, 3);
  scalar_mult(1/R, eRA, 3);
  scalar_mult(1/R, eRB, 3);

  /* compute polar and alpha angles */
  dot_arr(e12A, eRA, 3, &dot);
  theta_A = acos(dot);
  dot_arr(e12B, eRB, 3, &dot);
  theta_B = acos(dot);
  dot_arr(e32A, e12A, 3, &dot);
  alpha_A = acos(dot);
  dot_arr(e32B, e12B, 3, &dot);
  alpha_B = acos(dot);

  if (J==0) {      /* monomer-monomer distance or inverse distance */
    if (optinfo.frag_dist_rho)
      value = 1.0 / R / _bohr2angstroms;
    else
      value = R * _bohr2angstroms;
  }
  else if (J==1) { /* monomer polar angles */
    value = theta_A/_pi*180.0;
  }
  else if (J==2) {
    value = theta_B/_pi*180.0;
  }
  else if (J==3) { /* monomer-monomer torsion angle */
    cross_product(e12A,eRA,v);
    cross_product(e12B,eRA,v2);
    dot_arr(v, v2, 3, &dot);
    dot /= sin(theta_A) * sin(theta_B);
    value = acos(dot)/_pi*180.0;
    // determine sign
    cross_product(e12B,eRA,v);
    dot_arr(e12A, v, 3, &dot);
    if (dot < 0)
      value *= -1;
/* using the sin formula
    cross_product(e12B,eRA,v);
    dot_arr(e12A, v, 3, &dot);
    dot /= sin(theta_A) * sin(theta_B);
    value = asin(dot)/_pi*180.0;
*/
  }
  else if (J==4) { /* monomer-monomer internal rotation angles */
    cross_product(e32A,e12A,v);
    cross_product(e12A,eRA,v2);
    dot_arr(v, v2, 3, &dot);
    dot /= sin(theta_A) * sin(alpha_A);
    value = acos(dot)/_pi*180.0;
    // determine sign
    cross_product(eRA,e12A,v);
    dot_arr(e32A, v, 3, &dot);
    if (dot < 0)
      value *= -1;
/* using the sine formula
    cross_product(eRA,e12A,v);
    dot_arr(e32A, v, 3, &dot);
    dot /= sin(theta_A) * sin(alpha_A);
    value = asin(dot)/_pi*180.0;
*/
  }
  else if (J==5) {
    cross_product(e32B,e12B,v);
    cross_product(e12B,eRB,v2);
    dot_arr(v, v2, 3, &dot);
    dot /= sin(theta_B) * sin(alpha_B);
    value = acos(dot)/_pi*180.0;
    // determine sign
    cross_product(eRB,e12B,v);
    dot_arr(e32B, v, 3, &dot);
    if (dot < 0)
      value *= -1;
/* using the sine formula
    cross_product(eRB,e12B,v); // (e12B x eRA) = (eRB x e12B)
    dot_arr(e32B, v, 3, &dot);
    dot /= sin(theta_B) * sin(alpha_B);
    value = asin(dot)/_pi*180.0;
*/
  }

  free(e12A); free(e12B);
  free(e32A); free(e32B);
  free(eRA); free(eRB);
  free_block(dkA); free_block(dkB);
  free(v); free(v2); free(v3);
}

/** compute S vectors.  In this case, the derivative of the
    interfragment coordinates with respect to changes in the
    cartesian coordinates of the 6 reference atoms (3 on each fragment) **/

void fragment_class::compute_s(int natom, double *geom) {
  int xyz, k, i;
  double **dkA, **dkB; /* location of reference points */
  double *e12A, *e12B, *e32A, *e32B, *eRA, *eRB, *v, *v2, *geom_ang;
  double  d12A,  d12B,  d32A,  d32B,  R, c1, c2;
  double dot, alpha_A, alpha_B, theta_A, theta_B;

  geom_ang  = new double[3*natom];
  for (i=0;i<natom*3;++i)
    geom_ang[i] = geom[i] * _bohr2angstroms;

  dkA = block_matrix(3,3);
  dkB = block_matrix(3,3);
  e12A = init_array(3);
  e12B = init_array(3);
  e32A = init_array(3);
  e32B = init_array(3);
  eRA = init_array(3);
  eRB = init_array(3);
  v = init_array(3);
  v2 = init_array(3);

  /* compute reference points within each fragment */
  for (k=0; k<3; ++k) { /* k = point 1, 2 or 3 */
    for (xyz=0; xyz<3; ++xyz) {
      for (i=0; i<A_natom; ++i) 
        dkA[k][xyz] += A_weight[k*A_natom+i] * geom_ang[3*A_atom[i]+xyz];
      for (i=0; i<B_natom; ++i) 
        dkB[k][xyz] += B_weight[k*B_natom+i] * geom_ang[3*B_atom[i]+xyz];
    }
  }

  /* compute e vectors */
  for (xyz=0; xyz<3; ++xyz) {
    e12A[xyz] = dkA[1][xyz] - dkA[0][xyz];
    e12B[xyz] = dkB[1][xyz] - dkB[0][xyz];
    e32A[xyz] = dkA[1][xyz] - dkA[2][xyz];
    e32B[xyz] = dkB[1][xyz] - dkB[2][xyz];
    eRA[xyz]  = dkB[0][xyz] - dkA[0][xyz];
    eRB[xyz]  = dkA[0][xyz] - dkB[0][xyz];
  }
  d12A = sqrt( SQR(e12A[0]) +SQR(e12A[1]) +SQR(e12A[2]));
  d12B = sqrt( SQR(e12B[0]) +SQR(e12B[1]) +SQR(e12B[2]));
  d32A = sqrt( SQR(e32A[0]) +SQR(e32A[1]) +SQR(e32A[2]));
  d32B = sqrt( SQR(e32B[0]) +SQR(e32B[1]) +SQR(e32B[2]));
  R    = sqrt( SQR( eRA[0]) +SQR( eRA[1]) +SQR( eRA[2]));
  scalar_mult(1/d12A, e12A, 3);
  scalar_mult(1/d12B, e12B, 3);
  scalar_mult(1/d32A, e32A, 3);
  scalar_mult(1/d32B, e32B, 3);
  scalar_mult(1/R , eRA, 3);
  scalar_mult(1/R , eRB, 3);

  /* compute polar and alpha angles */
  dot_arr(e12A, eRA, 3, &dot);
  theta_A = acos(dot);
  dot_arr(e12B, eRA, 3, &dot);
  theta_B = acos(-1*dot);
  dot_arr(e32A, e12A, 3, &dot);
  alpha_A = acos(dot);
  dot_arr(e32B, e12B, 3, &dot);
  alpha_B = acos(dot);

  /* comments refer to counting atoms from 1 to match equations */
  if (J==0) { /* S vectors for coordinate 0 - distance - R or 1/R - atoms A1,B1 */
    for (xyz=0; xyz<3; ++xyz) {
      if (optinfo.frag_dist_rho) {
        A_s[3*0+xyz] = SQR(1.0/R) * eRA[xyz];
        B_s[3*0+xyz] = SQR(1.0/R) * eRB[xyz];
      }
      else {
        A_s[3*0+xyz] = -1.0 * eRA[xyz];
        B_s[3*0+xyz] = -1.0 * eRB[xyz];
      }
    }
  }
  else if (J==1) { /* S vectors for coordinate 1 - monomer polar angle on A */
    c1 = (d12A - R * cos(theta_A));
    c2 = (R - d12A * cos(theta_A));
    for (xyz=0; xyz<3; ++xyz) {
      A_s[3*0+xyz] = (c1 * e12A[xyz] + c2 * eRA[xyz]) / (d12A * R * sin(theta_A));
      A_s[3*1+xyz] = (cos(theta_A) * e12A[xyz] - eRA[xyz]) / (d12A * sin(theta_A));
      B_s[3*0+xyz] = (cos(theta_A) * eRA[xyz] - e12A[xyz]) / (R * sin(theta_A));
    }
  }
  else if (J==2) { /* S vectors for coordinate 2 - monomer polar angle on B */
    c1 = (d12B - R * cos(theta_B));
    c2 = (R - d12B * cos(theta_B));
    for (xyz=0; xyz<3; ++xyz) {
      B_s[3*0+xyz] = (c1 * e12B[xyz] + c2 * eRB[xyz]) / (d12B * R * sin(theta_B));
      B_s[3*1+xyz] = (cos(theta_B) * e12B[xyz] - eRB[xyz]) / (d12B * sin(theta_B));
      A_s[3*0+xyz] = (cos(theta_B) * eRB[xyz] - e12B[xyz]) / (R * sin(theta_B));
    }
  }
  else if (J==3) { /* S vectors for coordinate 3 - monomer-monomer torsion angle */
    cross_product(e12A, eRA, v);
    for (xyz=0; xyz<3; ++xyz)
      A_s[3*1+xyz] = v[xyz] / (d12A * SQR(sin(theta_A)));

    cross_product(e12B, eRB, v);
    for (xyz=0; xyz<3; ++xyz)
      B_s[3*1+xyz] = v[xyz] / (d12B * SQR(sin(theta_B)));
  
    cross_product(eRA, e12A, v);
    cross_product(eRA, e12B, v2);
    c1 = (R - d12A * cos(theta_A)) / (d12A * R * SQR(sin(theta_A)));
    c2 = cos(theta_B) / (R * SQR(sin(theta_B)));
    for (xyz=0; xyz<3; ++xyz)
      A_s[3*0+xyz] = c1 * v[xyz] - c2 * v2[xyz];
  
    cross_product(eRB, e12B, v);
    cross_product(eRB, e12A, v2);
    c1 = (R - d12B * cos(theta_B)) / (d12B * R * SQR(sin(theta_B)));
    c2 = cos(theta_A) / (R * SQR(sin(theta_A)));
    for (xyz=0; xyz<3; ++xyz)
      B_s[3*0+xyz] = c1 * v[xyz] - c2 * v2[xyz];
  }
  else if (J==4) { /* S vectors for coordinate 4 - internal rotation angle on A */
    /* atom 1 on A */
    cross_product(e12A, eRA, v);
    cross_product(e12A, e32A, v2);
    c1 = (d12A - R * cos(theta_A)) / (d12A * R * SQR(sin(theta_A)));
    c2 = cos(alpha_A) / (d12A * SQR(sin(alpha_A)));
    for (xyz=0; xyz<3; ++xyz)
      A_s[3*0+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 2 on A */
    cross_product(e12A, e32A, v);
    cross_product(e12A, eRA, v2);
    c1 = (d12A - d32A * cos(alpha_A)) / (d12A * d32A * SQR(sin(alpha_A)));
    c2 = cos(theta_A) / (d12A * SQR(sin(theta_A)));
    for (xyz=0; xyz<3; ++xyz)
      A_s[3*1+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 3 on A */
    cross_product(e32A, e12A, v);
    for (xyz=0; xyz<3; ++xyz)
      A_s[3*2+xyz] = v[xyz] / (d32A * SQR(sin(alpha_A)));

    /* atom 1 on B */
    cross_product(eRA, e12A, v);
    for (xyz=0; xyz<3; ++xyz)
      B_s[3*0+xyz] = v[xyz] / (R * SQR(sin(theta_A)));
  }
  else if (J==5) { /* S vectors for coordinate 5 - internal rotation angle on B */
    /* atom 1 on B */
    cross_product(e12B, eRB, v);
    cross_product(e12B, e32B, v2);
    c1 = (d12B - R * cos(theta_B)) / (d12B * R * SQR(sin(theta_B)));
    c2 = cos(alpha_B) / (d12B * SQR(sin(alpha_B)));
    for (xyz=0; xyz<3; ++xyz)
      B_s[3*0+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 2 on B */
    cross_product(e12B, e32B, v);
    cross_product(e12B, eRB, v2);
    c1 = (d12B - d32B * cos(alpha_B)) / (d12B * d32B * SQR(sin(alpha_B)));
    c2 = cos(theta_B) / (d12B * SQR(sin(theta_B)));
    for (xyz=0; xyz<3; ++xyz)
      B_s[3*1+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 3 on B */
    cross_product(e32B, e12B, v);
    for (xyz=0; xyz<3; ++xyz)
      B_s[3*2+xyz] = v[xyz] / (d32B * SQR(sin(alpha_B)));

    /* atom 1 on A */
    cross_product(eRB, e12B, v);
    for (xyz=0; xyz<3; ++xyz)
      A_s[3*0+xyz] = v[xyz] / (R * SQR(sin(theta_B)));
  }

  free(e12A); free(e12B);
  free(e32A); free(e32B);
  free(eRA); free(eRB);
  free(v); free(v2);
  free_block(dkA); free_block(dkB);
}

void fragment_class::print_s(void) {
  fprintf(outfile,"S_vectors on Fragment A ref points: DJ/dk\n");
  print_mat(&A_s, 1, 9, outfile);
  fprintf(outfile,"S_vectors on Fragment B ref points: DJ/dk\n");
  print_mat(&B_s, 1, 9, outfile);
}


