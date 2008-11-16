/*! \file
    \ingroup OPTKING
    \brief six coordinates for non-bonded fragments
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <physconst.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "fragment.h"

namespace psi { namespace optking {

void fragment_class::print(FILE *fp_out, int print_values, int print_weights) {
  int I,i,a,b;
  fprintf(fp_out,"\t(%d ", id);

  /* print_values implicitly lists which coordinates are present */
  if (!print_values) {
    fprintf(fp_out,"(");
    for (I=0; I<6; ++I)
      if (get_coord_on(I))
        fprintf(fp_out,"1 ");
      else
        fprintf(fp_out,"0 ");
    fprintf(fp_out,")");
  }

  fprintf(fp_out,"(", get_id());
  for (i=0; i<A_natom; ++i)
    fprintf(fp_out,"%d ", get_A_atom(i)+1);
  fprintf(outfile,")");
  fprintf(fp_out,"(", id);
  for (i=0; i<B_natom; ++i)
    fprintf(fp_out,"%d ", get_B_atom(i)+1);
  fprintf(outfile,")");

  if (print_weights) {
    fprintf(outfile,"\n\tFragment A reference point weights:\n");
    for (i=0; i<A_P; ++i) {
        fprintf(outfile,"\t(");
      for (a=0; a<A_natom; ++a)
        fprintf(outfile," %.2lf", get_A_weight(i,a));
      fprintf(outfile,")\n");
    }

    fprintf(outfile,"\tFragment B reference point weights:\n");
    for (i=0; i<B_P; ++i) {
        fprintf(outfile,"\t(");
      for (b=0; b<B_natom; ++b)
        fprintf(outfile," %.2lf", get_B_weight(i,b));
      fprintf(outfile,")\n");
    }
  }

  fprintf(outfile,"\t)\n");

  if (print_values) {
    if (coord_on[0]) {
      if (optinfo.frag_dist_rho)
        fprintf(fp_out, "\t\t 1/R(AB)  = %12.6lf\n", get_value(0));
      else
        fprintf(fp_out, "\t\t R(AB)    = %12.6lf\n", get_value(0));
    }
    if (coord_on[1])
      fprintf(fp_out,   "\t\t theta-A  = %12.6lf\n", get_value(1));
    if (coord_on[2])
      fprintf(fp_out,   "\t\t theta-B  = %12.6lf\n", get_value(2));
    if (coord_on[3])
      fprintf(fp_out,   "\t\t tau A-B  = %12.6lf\n", get_value(3));
    if (coord_on[4])
      fprintf(fp_out,   "\t\t chi-A    = %12.6lf\n", get_value(4));
    if (coord_on[5])
      fprintf(fp_out,   "\t\t chi-B    = %12.6lf\n", get_value(5));
  }
  else fprintf(outfile,"\n");
}

double fragment_class::get_val_A_or_rad(int I) {
  double tval;
  if (!coord_on[I])
    throw("fragment.get_val_A_or_rad() - coordinate is not active\n");
  if (I==0)
    tval = get_value(I);
  else if (I==1)
    tval = get_value(I) *_pi/180.0;
  else if (I==2)
    tval = get_value(I) *_pi/180.0;
  else if (I==3)
    tval = get_value(I) *_pi/180.0;
  else if (I==4)
    tval = get_value(I) *_pi/180.0;
  else if (I==5)
    tval = get_value(I) *_pi/180.0;
  return tval;
}

void fragment_class::compute(double *geom) {
  int xyz, k, i, a, I, b;
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
  for (k=0; k<A_P; ++k) { /* k = point 1, 2 or 3 */
    for (a=0; a<A_natom; ++a) {
      for (xyz=0; xyz<3; ++xyz) {
        dkA[k][xyz] += get_A_weight(k,a) * geom[3*get_A_atom(a)+xyz];
      }
    }
  }
  for (k=0; k<B_P; ++k) { /* k = point 1, 2 or 3 */
    for (b=0; b<B_natom; ++b) {
      for (xyz=0; xyz<3; ++xyz) {
        dkB[k][xyz] += get_B_weight(k,b) * geom[3*get_B_atom(b)+xyz];
      }
    }
  }

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

  if (coord_on[0]) {      /* monomer-monomer distance or inverse distance */
    if (optinfo.frag_dist_rho)
      value[0] = 1.0 / R / _bohr2angstroms;
    else
      value[0] = R * _bohr2angstroms;
  }
  if (coord_on[1]) { /* monomer polar angles */
    value[1] = theta_A/_pi*180.0;
  }
  if (coord_on[2]) {
    value[2] = theta_B/_pi*180.0;
  }
  if (coord_on[3]) { /* monomer-monomer torsion angle */
    cross_product(e12A,eRA,v);
    cross_product(e12B,eRA,v2);
    dot_arr(v, v2, 3, &dot);

    if ((sin(theta_A) > optinfo.sin_phi_denominator_tol) &&
        (sin(theta_B) > optinfo.sin_phi_denominator_tol)) {
         dot /= sin(theta_A);
         dot /= sin(theta_B);
    }
    else dot = 2.0 ;

//fprintf(outfile,"dot: %15.10lf\n",dot);

    if (dot > optinfo.cos_tors_near_1_tol) value[3] = 0.0 ;
    else if (dot < optinfo.cos_tors_near_neg1_tol) value[3] = 180.0 ;
    else {
//fprintf(outfile,"using acos to get value (J=3)\n");
      value[3] = acos(dot) / _pi * 180.0;
      // determine sign
      cross_product(e12B,eRA,v);
      dot_arr(e12A, v, 3, &dot);
      if (dot < 0) value[3] *= -1;
      /* using the sin formula
          cross_product(e12B,eRA,v);
          dot_arr(e12A, v, 3, &dot);
          dot /= sin(theta_A) * sin(theta_B);
          value[3] = asin(dot)/_pi*180.0; */
    }
  }
  if (coord_on[4]) { /* monomer-monomer internal rotation angles */
    cross_product(e32A,e12A,v);
    cross_product(e12A,eRA,v2);
    dot_arr(v, v2, 3, &dot);

    if ((sin(theta_A) > optinfo.sin_phi_denominator_tol) &&
        (sin(alpha_A) > optinfo.sin_phi_denominator_tol)) {
        dot /= sin(theta_A);
        dot /= sin(alpha_A);
    }
    else dot = 2.0;

    if (dot > optinfo.cos_tors_near_1_tol) value[4] = 0.0 ;
    else if (dot < optinfo.cos_tors_near_neg1_tol) value[4] = 180.0 ;
    else {
      value[4] = acos(dot) / _pi * 180.0;
      // determine sign
      cross_product(eRA,e12A,v);
      dot_arr(e32A, v, 3, &dot);
      if (dot < 0) value[4] *= -1;
      /* using the sine formula
          cross_product(eRA,e12A,v);
          dot_arr(e32A, v, 3, &dot);
          dot /= sin(theta_A) * sin(alpha_A);
          value[4] = asin(dot)/_pi*180.0; */
    }
  }
  if (coord_on[5]) {
    cross_product(e32B,e12B,v);
    cross_product(e12B,eRB,v2);
    dot_arr(v, v2, 3, &dot);
    if ((sin(theta_B) > optinfo.sin_phi_denominator_tol) &&
        (sin(alpha_B) > optinfo.sin_phi_denominator_tol)) {
      dot /= sin(theta_B);
      dot /= sin(alpha_B);
    }
    else
      dot = 2.0;

    if (dot > optinfo.cos_tors_near_1_tol) value[5] = 0.0 ;
    else if (dot < optinfo.cos_tors_near_neg1_tol) value[5] = 180.0 ;
    else {
      value[5] = acos(dot) / _pi * 180.0;
      // determine sign
      cross_product(eRB,e12B,v);
      dot_arr(e32B, v, 3, &dot);
      if (dot < 0) value[5] *= -1;
      /* using the sine formula
          cross_product(eRB,e12B,v); // (e12B x eRA) = (eRB x e12B)
          dot_arr(e32B, v, 3, &dot);
          dot /= sin(theta_B) * sin(alpha_B);
          value[5] = asin(dot)/_pi*180.0; */
    }
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
  int xyz, k, i, a, b;
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
  for (k=0; k<A_P; ++k)
    for (a=0; a<A_natom; ++a)
      for (xyz=0; xyz<3; ++xyz)
        dkA[k][xyz] += get_A_weight(k,a) * geom_ang[3*get_A_atom(a)+xyz];

  for (k=0; k<B_P; ++k)
    for (b=0; b<B_natom; ++b)
      for (xyz=0; xyz<3; ++xyz)
        dkB[k][xyz] += get_B_weight(k,b) * geom_ang[3*get_B_atom(b)+xyz];

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
  /* S vectors for coordinate 0 - distance - R or 1/R - atoms A1,B1 */
  if (coord_on[0]) {
    for (xyz=0; xyz<3; ++xyz) {
      if (optinfo.frag_dist_rho) {
        A_s[0][3*0+xyz] = SQR(1.0/R) * eRA[xyz];
        B_s[0][3*0+xyz] = SQR(1.0/R) * eRB[xyz];
      }
      else {
        A_s[0][3*0+xyz] = -1.0 * eRA[xyz];
        B_s[0][3*0+xyz] = -1.0 * eRB[xyz];
      }
    }
  }

  /* S vectors for coordinate 1 - polar angle A2, A1, B1 */
  if (coord_on[1]) {
    c1 = (d12A - R * cos(theta_A));
    c2 = (R - d12A * cos(theta_A));
    for (xyz=0; xyz<3; ++xyz) {
      A_s[1][3*0+xyz] = (c1 * e12A[xyz] + c2 * eRA[xyz]) / (d12A * R * sin(theta_A));
      A_s[1][3*1+xyz] = (cos(theta_A) * e12A[xyz] - eRA[xyz]) / (d12A * sin(theta_A));
      B_s[1][3*0+xyz] = (cos(theta_A) * eRA[xyz] - e12A[xyz]) / (R * sin(theta_A));
    }
  }

  /* S vectors for coordinate 2 - polar angle B2, B1, A1 */
  if (coord_on[2]) {
    c1 = (d12B - R * cos(theta_B));
    c2 = (R - d12B * cos(theta_B));
    for (xyz=0; xyz<3; ++xyz) {
      B_s[2][3*0+xyz] = (c1 * e12B[xyz] + c2 * eRB[xyz]) / (d12B * R * sin(theta_B));
      B_s[2][3*1+xyz] = (cos(theta_B) * e12B[xyz] - eRB[xyz]) / (d12B * sin(theta_B));
      A_s[2][3*0+xyz] = (cos(theta_B) * eRB[xyz] - e12B[xyz]) / (R * sin(theta_B));
    }
  }

  /* S vectors for coordinate 3 - tau torsion B2, B1, A1 A2 */
  if (coord_on[3]) {
    cross_product(e12A, eRA, v);
    for (xyz=0; xyz<3; ++xyz)
      A_s[3][3*1+xyz] = v[xyz] / (d12A * SQR(sin(theta_A)));

    cross_product(e12B, eRB, v);
    for (xyz=0; xyz<3; ++xyz)
      B_s[3][3*1+xyz] = v[xyz] / (d12B * SQR(sin(theta_B)));
  
    cross_product(eRA, e12A, v);
    cross_product(eRA, e12B, v2);
    c1 = (R - d12A * cos(theta_A)) / (d12A * R * SQR(sin(theta_A)));
    c2 = cos(theta_B) / (R * SQR(sin(theta_B)));
    for (xyz=0; xyz<3; ++xyz)
      A_s[3][3*0+xyz] = c1 * v[xyz] - c2 * v2[xyz];
  
    cross_product(eRB, e12B, v);
    cross_product(eRB, e12A, v2);
    c1 = (R - d12B * cos(theta_B)) / (d12B * R * SQR(sin(theta_B)));
    c2 = cos(theta_A) / (R * SQR(sin(theta_A)));
    for (xyz=0; xyz<3; ++xyz)
      B_s[3][3*0+xyz] = c1 * v[xyz] - c2 * v2[xyz];
  }

  /* S vectors for coordinate 4 - chi_A: A3, A2, A1, B1 */
  if (coord_on[4]) {
    /* atom 1 on A */
    cross_product(e12A, eRA, v);
    cross_product(e12A, e32A, v2);
    c1 = (d12A - R * cos(theta_A)) / (d12A * R * SQR(sin(theta_A)));
    c2 = cos(alpha_A) / (d12A * SQR(sin(alpha_A)));
    for (xyz=0; xyz<3; ++xyz)
      A_s[4][3*0+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 2 on A */
    cross_product(e12A, e32A, v);
    cross_product(e12A, eRA, v2);
    c1 = (d12A - d32A * cos(alpha_A)) / (d12A * d32A * SQR(sin(alpha_A)));
    c2 = cos(theta_A) / (d12A * SQR(sin(theta_A)));
    for (xyz=0; xyz<3; ++xyz)
      A_s[4][3*1+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 3 on A */
    cross_product(e32A, e12A, v);
    for (xyz=0; xyz<3; ++xyz)
      A_s[4][3*2+xyz] = v[xyz] / (d32A * SQR(sin(alpha_A)));

    /* atom 1 on B */
    cross_product(eRA, e12A, v);
    for (xyz=0; xyz<3; ++xyz)
      B_s[4][3*0+xyz] = v[xyz] / (R * SQR(sin(theta_A)));
  }

  /* S vectors for coordinate 5 - chi_B: B3, B2, B1, A1 */
  if (coord_on[5]) {
    /* atom 1 on B */
    cross_product(e12B, eRB, v);
    cross_product(e12B, e32B, v2);
    c1 = (d12B - R * cos(theta_B)) / (d12B * R * SQR(sin(theta_B)));
    c2 = cos(alpha_B) / (d12B * SQR(sin(alpha_B)));
    for (xyz=0; xyz<3; ++xyz)
      B_s[5][3*0+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 2 on B */
    cross_product(e12B, e32B, v);
    cross_product(e12B, eRB, v2);
    c1 = (d12B - d32B * cos(alpha_B)) / (d12B * d32B * SQR(sin(alpha_B)));
    c2 = cos(theta_B) / (d12B * SQR(sin(theta_B)));
    for (xyz=0; xyz<3; ++xyz)
      B_s[5][3*1+xyz] = c1 * v[xyz] + c2 * v2[xyz];

    /* atom 3 on B */
    cross_product(e32B, e12B, v);
    for (xyz=0; xyz<3; ++xyz)
      B_s[5][3*2+xyz] = v[xyz] / (d32B * SQR(sin(alpha_B)));

    /* atom 1 on A */
    cross_product(eRB, e12B, v);
    for (xyz=0; xyz<3; ++xyz)
      A_s[5][3*0+xyz] = v[xyz] / (R * SQR(sin(theta_B)));
  }

  free(e12A); free(e12B);
  free(e32A); free(e32B);
  free(eRA); free(eRB);
  free(v); free(v2);
  free_block(dkA); free_block(dkB);
}

void fragment_class::print_s(void) {
  int i,I;
  fprintf(outfile,"S_vectors on fragment reference points: DJ/dk\n");
  for (I=0; I<6; ++I) {
    if (coord_on[I]) {
      fprintf(outfile,"A:");
      for (i=0;i<(3*A_P);++i)
        fprintf(outfile," %11.7lf", A_s[I][i]);
      fprintf(outfile,"\n");
      fprintf(outfile,"B:");
      for (i=0;i<(3*B_P);++i)
        fprintf(outfile," %11.7lf", B_s[I][i]);
      fprintf(outfile,"\n");
    }
  }
}

}} /* namespace psi::optking */
