#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <file30.h>
#include <qt.h>
#include <math.h>
#include "globaldefs.h"
#include "globals.h"

#define MO_HESS_MIN 1.0E-2


/*
** calc_orb_step()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and an approximate orbital Hessian
**
** C. David Sherrill
** April 1998
*/
void calc_orb_step(int npairs, double *grad, double *hess_diag, double *theta)
{

  int pair;
  double numer, denom;

  for (pair=0; pair<npairs; pair++) {
    numer = grad[pair];
    denom = hess_diag[pair];
    if (denom < 0.0) {
      fprintf(outfile, "Warning: MO Hessian denominator negative\n");
      denom = -denom;
    }
    if (denom < MO_HESS_MIN) {
      fprintf(outfile, "Warning: MO Hessian denominator too small\n");
      denom = MO_HESS_MIN;
    } 
    theta[pair] =  - numer / denom;
  }

}


/*
** print_step
**
** This function prints out the information for a given orbital iteration
*/
int print_step(int npairs, int steptype)
{
  FILE *sumfile;
  char sumfile_name[] = "file14.dat";
  int i, entries, iter, *nind;
  double *rmsgrad, *scaled_rmsgrad, *energies, energy;
  char **comments;

  /* open ascii file, get number of entries already in it */

  sumfile = fopen(sumfile_name, "r");
  if (sumfile == NULL) { /* the file doesn't exist yet */
    entries = 0;
    if (Params.print_lvl)
      fprintf(outfile, "\nPreparing new file %s\n", sumfile_name);
  }
  else {
    if (fscanf(sumfile, "%d", &entries) != 1) {
      fprintf(outfile,"(print_step): Trouble reading num entries in file %s\n",
        sumfile_name);
      fclose(sumfile);
      return;
    }
  }

  rmsgrad = init_array(entries+1);
  scaled_rmsgrad = init_array(entries+1);
  energies= init_array(entries+1);
  nind = init_int_array(entries+1);
  comments = (char **) malloc ((entries+1) * sizeof (char *));
  for (i=0; i<entries+1; i++) {
    comments[i] = (char *) malloc (MAX_COMMENT * sizeof(char));
  }

  for (i=0; i<entries; i++) {
    fscanf(sumfile, "%d %d %lf %lf %lf %s", &iter, &(nind[i]), 
           &(scaled_rmsgrad[i]), &(rmsgrad[i]), &(energies[i]), comments[i]);
  }

  file30_init();
  energy = file30_rd_ecorr();
  file30_close();

  scaled_rmsgrad[entries] = CalcInfo.scaled_mo_grad_rms;
  rmsgrad[entries] = CalcInfo.mo_grad_rms;
  energies[entries] = energy;
  nind[entries] = npairs;

  if (steptype == 0) 
    strcpy(comments[entries], "CONV");
  else if (steptype == 1)
    strcpy(comments[entries], "NR");
  else if (steptype == 2)
    strcpy(comments[entries], "DIIS"); 
  else {
    fprintf(outfile, "(print_step): Unrecognized steptype %d\n", steptype);
    strcpy(comments[entries], "?");
  }

  if (entries) fclose(sumfile);

  /* now open file for writing, write out old info plus new */
  if ((sumfile = fopen("file14.dat", "w")) == NULL) {
    fprintf(outfile, "(print_step): Unable to open file %s\n", sumfile_name);
  }
  else {
    entries++;
    fprintf(sumfile, "%5d\n", entries);
    for (i=0; i<entries; i++) {
      fprintf(sumfile, "%5d %5d %14.9lf %14.9lf %20.12lf %9s\n", i+1, nind[i], 
              scaled_rmsgrad[i], rmsgrad[i], energies[i], comments[i]);
    }
    fclose(sumfile);
  }

  free(scaled_rmsgrad);
  free(rmsgrad);
  free(energies);
  free(nind);
  for (i=0; i<entries; i++)
    free(comments[i]);
  free(comments);

}

