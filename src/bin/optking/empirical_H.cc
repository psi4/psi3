// This function generates an empirical guess Hessian from a given set of
//   salcs according to Schlegel, Theor. Chim. Acta, 66, 333 (1984).

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
  #include <stdio.h>
  #include <libchkpt/chkpt.h>
  #include <stdlib.h>
  #include <string.h>
  #include <libciomr/libciomr.h>
  #include <physconst.h>
  #include <libipv1/ip_lib.h>
  #include <physconst.h>
  #include <libpsio/psio.h>
  #include <psifiles.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

static double cov_radii[37] = {
       0.32, 0.32, 0.60,
       1.2, 1.05,        0.81, 0.77, 0.74, 0.74, 0.72, 0.72,
       1.5, 1.4,         1.3,  1.17, 1.10, 1.04, 0.99, 0.99,
       1.8, 1.6,
          1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4,
                         1.4,  1.3,  1.2,  1.2,  1.1,  1.1};

void empirical_H(internals &simples, salc_set &symm, cartesians &carts) {
   int i, j, k, atomA, atomB, atomC, atomD, simple, count = 0;
   int perA, perB, col;
   double A,B,C, rBC, r1[3], r2[3], r3[3];
   double *f, val, cov_radius, *coord, norm_r1, norm_r2, norm_r3;

   f = init_array(simples.get_num());      
   coord = carts.get_coord();

// Form diagonal Hessian in simple internals first
   count = -1;
   for (i=0;i<simples.stre.get_num();++i) {
      atomA = simples.stre.get_A(i);
      atomB = simples.stre.get_B(i);
      A = 1.734;
      B = 1.0;
      perA = 1;
      perB = 1;
      /*fprintf(outfile,"\n\n stretch given %d and %d \n",atomA,atomB);*/
      if ((int) (carts.get_atomic_num(atomA)) == 1)
         perA = 1;
      else if ((int) (carts.get_atomic_num(atomA)) <= 10)
         perA = 2;
      else if ((int) (carts.get_atomic_num(atomA)) <= 18)
         perA = 3;
      else if ((int) (carts.get_atomic_num(atomA)) <= 36)
         perA = 4;
      else
         perA = 5;
    
      if ((int) (carts.get_atomic_num(atomB)) == 1)
         perB = 1;
      else if ((int) (carts.get_atomic_num(atomB)) <= 10)
         perB = 2;
      else if ((int) (carts.get_atomic_num(atomB)) <= 18)
         perB = 3;
      else if ((int) (carts.get_atomic_num(atomB)) <= 36)
         perB = 4;
      else
         perB = 5;
   
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
      val = simples.stre.get_val(i) / _bohr2angstroms;
      val = A/((val-B)*(val-B)*(val-B));
      f[++count] = val * _hartree2J*1.0E18/SQR(_bohr2angstroms);
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
      atomA = simples.tors.get_A(i);
      atomB = simples.tors.get_B(i);
      atomC = simples.tors.get_C(i);
      atomD = simples.tors.get_D(i);
      A = 0.0023;
      B = 0.07;
      if (  ((int) (carts.get_atomic_num(atomB)) < 37)
        &&  ((int) (carts.get_atomic_num(atomC)) < 37) )
        cov_radius = (cov_radii[(int) (carts.get_atomic_num(atomB))]
         + cov_radii[(int) (carts.get_atomic_num(atomC))])/_bohr2angstroms;
      else
        cov_radius = 1.1/_bohr2angstroms;

      rBC = 0.0;
      for (j=0;j<3;++j)
         rBC += SQR(coord[3*atomB+j] - coord[3*atomC+j]); 
      rBC = sqrt(rBC);
      if (rBC > (cov_radius + A/B))
         B = 0.0;
      f[++count] = _hartree2J*1.0E18 * (A - (B*(rBC - cov_radius)));
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
      atomA = simples.lin_bend.get_A(i);
      atomB = simples.lin_bend.get_B(i);
      atomC = simples.lin_bend.get_C(i);
      val = 0.10;
      f[++count] = val * _hartree2J*1.0E18;
   }
   free(coord);

  // Now transform into delocalized coordinates U^t H U
   double **intcos;
   double **f_new;

   intcos = block_matrix(symm.get_num(),simples.get_num());
   int id, index;

   for (i=0;i<symm.get_num();++i) {
     for (j=0;j<symm.get_length(i);++j) {
       id = symm.get_simple(i,j);
       index = simples.id_to_index(id);
       intcos[i][index] = symm.get_coeff(i,j);
     }
   }

   f_new = block_matrix(symm.get_num(),symm.get_num());
   for (i=0;i<symm.get_num();++i)
     for (j=0;j<symm.get_num();++j)
       for (k=0;k<simples.get_num();++k)
         f_new[i][j] += intcos[i][k] * f[k] * intcos[j][k]; 

/*
   fprintf(outfile,"\nempirical force constants made\n");
   for (i=0;i<symm.get_num();++i) {
     for (j=0;j<symm.get_num();++j)
        fprintf(outfile,"%15.8lf",f_new[i][j]);
     fprintf(outfile,"\n");
   }
*/


  /*** write to PSIF_OPTKING ***/
   open_PSIF();
   psio_write_entry(PSIF_OPTKING, "Force Constants",
       (char *) &(f_new[0][0]),symm.get_num()*symm.get_num()*sizeof(double));
   close_PSIF();

/*
  // Print out forces constants to fconst.dat 
   fp_fconst = fopen("fconst.dat","w");
   for (i=0;i<symm.get_num();++i) {
     col = 0;
     for (j=0; j<=i ; ++j) {
       if (col == 8) {
         fprintf(fp_fconst,"\n");
         col = 0;
       }
       fprintf(fp_fconst,"%10.6f",f_new[i][j]);
       ++col;
     }
     fprintf(fp_fconst,"\n");
   }
   fclose(fp_fconst);
*/
   free(f);
   free_block(f_new);
   free_block(intcos);
   return;
}

