/*
** CALC_DE_DT
**
** This file contains the code to evaluate dE/dTheta  = dE/dU * dU/dTheta
**
** C. David Sherrill
** May 1998
*/

#include <stdlib.h>
#include <stdio.h>


/*
** calc_dE_dT()
**
** This function gives us dE/dTheta from dE/dU * dU/dTheta
**
** Feed in a vector of thetas (theta) and a matrix dE/dU and write
**   to a matrix of dE/d(theta)  (dET)
**
** Based on similar code from an old version of the VBD code.
** Can use blas calls to drot later.
**
*/
void calc_dE_dT(double **theta, double **dEU, double *dET)
{
  int i,a,a2,m,l;
  int nocc, nvir, n, num;
  double costheta, sintheta, temp;
  double **Uleft, **Uright, **Scratch;

  nocc = CalcInfo.nocc;  /* number occupied */
  nvir = CalcInfo.nvir;  /* number virtual */
  n = CalcInfo.nspcorr;  /* total number spin orbitals */

  /*
  fprintf(outfile, "Theta matrix\n");
  print_mat(theta, nocc, nvir, outfile);
  */

  zero_arr(dET, nocc*nvir);


  /* form temporary U matrices */

  Uleft = block_matrix(n ,n);
  Uright = block_matrix(n, n);
  Scratch = block_matrix(n, n);

  /* init Uright and Uleft to unit matrix */
  for(m=0; m < n; m++)  {
    Uright[m][m] = 1.0;
    Uleft[m][m] = 1.0;
  }

  /* init Uleft to the U matrix */
  for(i=0; i < nocc; i++)  {
    for(a=nocc,a2=0; a < n; a++,a2++)  {
      costheta = cos(theta[i][a2]);
      sintheta = sin(theta[i][a2]);

      /* in-place postmultication of Uleft by G_{ia} */
      for(m=0; m < n; m++)  {
        temp = Uleft[m][i];
        Uleft[m][i] = temp*costheta - Uleft[m][a]*sintheta;
        Uleft[m][a] = Uleft[m][a]*costheta + temp*sintheta;
      }
    }
  }
  /* give old dE/dU and returns 1/4
  **   backtransformed trnasposed dE/dU */

  /*
  fprintf(outfile, "dE/dU before backtransform: \n");
  print_mat(dEU, n, n, outfile);
  fprintf(outfile, "\nU matrix for transform: \n");
  print_mat(Uleft, n, n, outfile);
  fflush(outfile);

  */
  /* test out uleft */
  /* form_u(Uleft, nocc, nvir, n, theta); */

  calc_back_trans(dEU, Uleft);

  /*
  fprintf(outfile, "dE/dU after backtransform: \n");
  print_mat(dEU, n, n, outfile);
  fprintf(outfile, "Now entering derivative routine: \n");
  fflush(outfile);
  */

  /* Loop over i,a to form dE/d(theta): we are working right to left in this
     algorithm, hence we must go backwards through the theta list */
  for(i=nocc-1,num=nocc*nvir-1; i >= 0; i--) {
    for(a=n-1,a2=nvir-1; a >= nocc; a--,num--,a2--) {


      costheta = cos(theta[i][a2]);
      sintheta = sin(theta[i][a2]);

      /*
      fprintf(outfile, "Theta matrix\n");
      print_mat(theta, nocc, nvir, outfile);
      fprintf(outfile, "Derivative (i=%d, a=%d)\n", i, a);
      fprintf(outfile, "Cos = %lf, Sin=%lf\n", costheta, sintheta);
      */

      /* post-multiply Uleft by G(+) */
      for(m=0; m < n; m++)  {
        temp = Uleft[m][i];
        Uleft[m][i] = temp*costheta + Uleft[m][a]*sintheta;
        Uleft[m][a] = Uleft[m][a]*costheta - temp*sintheta;
      }

      /*
      fprintf(outfile, "Uleft after postmultiplication by G(+)(%d,%d)\n",
        i, a);
      print_mat(Uleft, n, n, outfile);
      */

      /* Now do Uleft*dG/d(theta)*Uright series of multi */
      for(l=0; l < n; l++)  {
        for(m=0; m < n; m++)  {
          Scratch[l][m] = (-sintheta*Uleft[l][i] - costheta*Uleft[l][a])
            *Uright[i][m] + (costheta*Uleft[l][i] - sintheta*Uleft[l][a])
              *Uright[a][m];
        }
      }

      /*
      fprintf(outfile, "Uleft * dG/dTheta(%d,%d) * Uright\n", i, a);
      print_mat(Scratch, n, n, outfile);
      */

      for(l=0; l < n; l++)  {
        for(m=0; m < n; m++)  {
          dET[num] += dEU[l][m]*Scratch[l][m];
          /*
          if ((fabs(Scratch[l][m]) > 0.0001) && ((l<nocc && m<nocc) ||
            (l >= nocc && m >= nocc))) {
            fprintf(outfile, "nonzero element for theta(%d, %d) element", i, a);
            fprintf(outfile, "%d %d\n", l, m);
          }
          */
        }
      }

      /*
      fprintf(outfile, "dE/dTheta(%d,%d) = %12.6lf\n", i, a, dET[num]);
      */

      /* pre-multiply Uright by G */
      for(m=0; m < n; m++)  {
        temp = Uright[i][m];
        Uright[i][m] = temp*costheta + Uright[a][m]*sintheta;
        Uright[a][m] = Uright[a][m]*costheta - temp*sintheta;
      }

      /*
      fprintf(outfile, "Uright after premultiplication by G \n");
      print_mat(Uright, n, n, outfile);
      */

    }
  }
  /* free memory */

  free_block(Uleft);
  free_block(Uright);
  free_block(Scratch);
}



