/*! \file
    \ingroup OPTKING
    \brief generates an empirical guess Hessian from a given set of
      salcs according to Schlegel, Theor. Chim. Acta, 66, 333 (1984) or
      according to Fischer and Almlof, J. Phys. Chem., 96, 9770 (1992).
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <physconst.h>
#include <cov_radii.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

namespace psi { namespace optking {

// returns covalent bond length in bohr from atomic numbers
inline double Rcov(int ZA, int ZB) {
  return (cov_radii[ZA] + cov_radii[ZB]) / _bohr2angstroms;
} 
inline double Rcov(double ZA, double ZB) {
  return (cov_radii[(int) ZA] + cov_radii[(int) ZB]) / _bohr2angstroms;
}

// returns period from atomic number
inline int period(int ZA) {
  if      (ZA ==  1) return 1;
  else if (ZA <= 10) return 2;
  else if (ZA <= 18) return 3;
  else if (ZA <= 36) return 4;
  else               return 5;
}
inline int period(double ZA) {
  if      ((int) ZA ==  1) return 1;
  else if ((int) ZA <= 10) return 2;
  else if ((int) ZA <= 18) return 3;
  else if ((int) ZA <= 36) return 4;
  else                     return 5;
}

// compute force constants between fragments
double fragment_distance_fc(fragment_class &frag, cartesians &carts);

void empirical_H(internals &simples, salc_set &symm, cartesians &carts) {
  int i, j, k, atomA, atomB, atomC, atomD, simple, count = -1, perA, perB, L;
  int a, b;
  double rAB, rBC, rABcov, rBCcov, rBDcov, prefactor_i;
  double A, B, C, D, E, r1[3], r2[3], r3[3], tval;
  double *f, val, *coord, norm_r1, norm_r2, norm_r3;

  f = init_array(simples.get_num());      
  coord = carts.get_coord();

  // Form diagonal Hessian in simple internals first
  if (optinfo.empirical_H == OPTInfo::SCHLEGEL) {
    for (i=0;i<simples.stre.get_num();++i) {
      atomA = simples.stre.get_A(i);
      atomB = simples.stre.get_B(i);
      perA = period(carts.get_atomic_num(atomA));
      perB = period(carts.get_atomic_num(atomB));
      A = 1.734;

      if ((perA==1) && (perB==1))
         B = -0.244;
      else if (((perA==1) && (perB==2)) || ((perB==1) && (perA==2)))
         B = 0.352;
      else if ((perA==2) && (perB==2))
         B = 1.085;
      else if (((perA==1) && (perB==3)) || ((perB==1) && (perA==3)))
         B = 0.660;
      else if (((perA==2) && (perB==3)) || ((perB==2) && (perA==3)))
         B = 1.522;
      else
         B = 2.068;

      rAB = carts.R(atomA,atomB);
      // fc in au / bohr^2
      tval = A/((rAB-B)*(rAB-B)*(rAB-B));
      // fc in aJ/Ang^2
      f[++count] = tval * _hartree2J*1.0E18/SQR(_bohr2angstroms);
    }
    for (i=0;i<simples.bend.get_num();++i) {
      atomA = simples.bend.get_A(i);
      atomB = simples.bend.get_B(i);
      atomC = simples.bend.get_C(i);
      if ( ((int) (carts.get_atomic_num(atomA)) == 1) ||
           ((int) (carts.get_atomic_num(atomC) == 1)) )
        val = 0.160;
      else
        val = 0.250;
      f[++count] = val * _hartree2J*1.0E18;
    }
    for (i=0;i<simples.tors.get_num();++i) {
      atomB = simples.tors.get_B(i);
      atomC = simples.tors.get_C(i);
      A = 0.0023;
      B = 0.07;
      rBCcov = Rcov(carts.get_atomic_num(atomB), carts.get_atomic_num(atomC));
      rBC = carts.R(atomB,atomC);
      if (rBC > (rBCcov + A/B)) B = 0.0; // keep > 0
      f[++count] = (A - (B*(rBC - rBCcov))) * _hartree2J*1.0E18;
    }
    for (i=0;i<simples.out.get_num();++i) {
      atomA = simples.out.get_A(i);
      atomB = simples.out.get_B(i);
      atomC = simples.out.get_C(i);
      atomD = simples.out.get_D(i);
      A = 0.045;
      for (j=0;j<3;++j) {
         r1[j] = coord[3*atomA+j] - coord[3*atomB+j];
         r2[j] = coord[3*atomC+j] - coord[3*atomB+j];
         r3[j] = coord[3*atomD+j] - coord[3*atomB+j];
      }
      norm_r1 = norm_r2 = norm_r3 = 0.0;
      for (j=0;j<3;++j) {
         norm_r1 += SQR(coord[3*atomA+j] - coord[3*atomB+j]);
         norm_r2 += SQR(coord[3*atomC+j] - coord[3*atomB+j]);
         norm_r3 += SQR(coord[3*atomD+j] - coord[3*atomB+j]);
      }
      norm_r1 = sqrt(norm_r1);
      norm_r2 = sqrt(norm_r2);
      norm_r3 = sqrt(norm_r3);
      val =  ( 1 - ((  r1[0]*r2[1]*r3[2]
                      +r1[2]*r2[0]*r3[1]
                      +r1[1]*r2[2]*r3[0]
                      -r1[2]*r2[1]*r3[0]
                      -r1[0]*r2[2]*r3[1]
                      -r1[1]*r2[0]*r3[2] ) / (norm_r1*norm_r2*norm_r3)));
      f[++count] = _hartree2J*1.0E18 * (A * val * val * val * val);
    }
    for (i=0;i<simples.lin_bend.get_num();++i) {
      f[++count] = 0.10 * _hartree2J*1.0E18;
    }
    for (i=0;i<simples.frag.get_num();++i) {
      if (simples.frag.get_coord_on(i,0)) {
        // if (optinfo.frag_dist_rho == 1) val = 100.0;
        // else val = 0.01;
        int min_a,min_b;
        rAB = 1e6;
        // find minimum distance between two atoms, one in each fragment
        for (a=0; a<simples.frag.get_A_natom(i); ++a) {
          atomA = simples.frag.get_A_atom(i,a);
          for (b=0; b<simples.frag.get_B_natom(i); ++b) {
            atomB = simples.frag.get_B_atom(i,b);
            tval = carts.R(atomA,atomB);
            if (tval < rAB) {
              rAB = tval;
              min_a = atomA;
              min_b = atomB;
            }
          }
        }
        rABcov = Rcov(carts.get_atomic_num(min_a), carts.get_atomic_num(min_b));
        A = 0.3601; B = 1.944;
        tval = A * exp(-B*(rAB - rABcov));
        if (optinfo.frag_dist_rho == 1)
          tval *= pow(rAB,4);
        tval *= _hartree2J*1.0E18 / SQR(_bohr2angstroms);
        f[++count] = tval;
      }
      if (simples.frag.get_coord_on(i,1))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,2))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,3))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,4))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,5))
        f[++count] = 0.001;
    }
  }
  else if (optinfo.empirical_H == OPTInfo::FISCHER) {
    for (i=0;i<simples.stre.get_num();++i) {
      atomA = simples.stre.get_A(i);
      atomB = simples.stre.get_B(i);

      rAB   = carts.R(atomA,atomB);
      rABcov = Rcov(carts.get_atomic_num(atomA), carts.get_atomic_num(atomB));

      A = 0.3601; B = 1.944;

      tval = A * exp(-B*(rAB - rABcov));
      f[++count] = tval * _hartree2J*1.0E18 / SQR(_bohr2angstroms);
    }
    for (i=0;i<simples.bend.get_num();++i) {
      atomA = simples.bend.get_A(i);
      atomB = simples.bend.get_B(i);
      atomC = simples.bend.get_C(i);

      rABcov = Rcov(carts.get_atomic_num(atomA), carts.get_atomic_num(atomB));
      rBCcov = Rcov(carts.get_atomic_num(atomB), carts.get_atomic_num(atomC));
      rAB = carts.R(atomA, atomB);
      rBC = carts.R(atomB, atomC);

      A = 0.089; B = 0.11; C = 0.44; D = -0.42;

      tval = A + B/pow(rABcov*rBCcov, D) * exp(-C*( rAB + rBC - rABcov - rBCcov));
      f[++count] = tval * _hartree2J*1.0E18;
    }
    for (i=0;i<simples.tors.get_num();++i) {
      atomB = simples.tors.get_B(i);
      atomC = simples.tors.get_C(i);

      rBCcov = Rcov(carts.get_atomic_num(atomB), carts.get_atomic_num(atomC));
      rBC = carts.R(atomB,atomC);

      A = 0.0015; B = 14.0; C = 2.85; D = 0.57; E = 4.00;
      L = 2; // # of bonds connected to atom a and b except the central bond?

      tval = B * pow(L,D) / pow(rBC * rBCcov, E);
      tval = A + tval * exp(-C * (rBC - rBCcov));
      f[++count] = tval * _hartree2J*1.0e18;
    }
    for (i=0;i<simples.out.get_num();++i) {
      atomA = simples.out.get_A(i);
      atomB = simples.out.get_B(i);
      atomC = simples.out.get_C(i);
      atomD = simples.out.get_D(i);
      val = simples.out.get_val(i)*_pi/180.0;

      A = 0.0025; B = 0.0061; C = 3.00; D = 4.00; E = 0.80;

      rABcov = Rcov(carts.get_atomic_num(atomA), carts.get_atomic_num(atomB));
      rBCcov = Rcov(carts.get_atomic_num(atomB), carts.get_atomic_num(atomC));
      rBDcov = Rcov(carts.get_atomic_num(atomB), carts.get_atomic_num(atomD));
      carts.R(atomA, atomB);

      tval = B * pow(rBCcov*rBDcov, E) * pow(cos(val),D);
      tval = A + tval * exp(-C*(rAB - rABcov));
      f[++count] = tval * _hartree2J*1.0E18 ;
    }
    for (i=0;i<simples.lin_bend.get_num();++i) {
      f[++count] = 0.10 * _hartree2J*1.0E18;
    }
    for (i=0;i<simples.frag.get_num();++i) {
      if (simples.frag.get_coord_on(i,0)) {
        int min_a,min_b;
        rAB = 1e6;
        // find minimum distance between two atoms, one in each fragment
        for (a=0; a<simples.frag.get_A_natom(i); ++a) {
          atomA = simples.frag.get_A_atom(i,a);
          for (b=0; b<simples.frag.get_B_natom(i); ++b) {
            atomB = simples.frag.get_B_atom(i,b);
            tval = carts.R(atomA,atomB);
            if (tval < rAB) {
              rAB = tval;
              min_a = atomA;
              min_b = atomB;
            }
          }
        }
        rABcov = Rcov(carts.get_atomic_num(min_a), carts.get_atomic_num(min_b));
        A = 0.3601; B = 1.944;
        tval = A * exp(-B*(rAB - rABcov));
        if (optinfo.frag_dist_rho == 1)
          tval *= pow(rAB,4);
        tval *= _hartree2J*1.0E18 / SQR(_bohr2angstroms);
        f[++count] = tval;
      }
      if (simples.frag.get_coord_on(i,1))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,2))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,3))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,4))
        f[++count] = 0.001;
      if (simples.frag.get_coord_on(i,5))
        f[++count] = 0.001;
    }
  }
  free(coord);

  //fprintf(outfile,"Diagonal force constants for simple internals\n");
  //print_mat(&f,1,simples.get_num(),outfile);

 // Now transform into salc coordinates U^t H U
  double **intcos;
  double **f_new;

  intcos = block_matrix(symm.get_num(),simples.get_num());
  int id, index;

  for (i=0;i<symm.get_num();++i) {
    prefactor_i = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) {
      id = symm.get_simple(i,j);
      index = simples.id_to_index(id);
      intcos[i][index] = prefactor_i * symm.get_coeff(i,j);
    }
  }

  // fprintf(outfile,"Simples to Salc matrix\n");
  // print_mat(intcos,symm.get_num(),simples.get_num(),outfile);

  f_new = block_matrix(symm.get_num(),symm.get_num());
  for (i=0;i<symm.get_num();++i)
    for (j=0;j<symm.get_num();++j)
      for (k=0;k<simples.get_num();++k)
        f_new[i][j] += intcos[i][k] * f[k] * intcos[j][k]; 

 /*** write to PSIF_OPTKING ***/
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(f_new[0][0]),symm.get_num()*symm.get_num()*sizeof(double));
  close_PSIF();

  free(f);
  free_block(f_new);
  free_block(intcos);
  return;
}

}}
