/*
**
** DETCASMAN
**
** Program to manage the iteration of (transqt, detci, clag, detcas)
** required for orbital optimization using the DETCAS program.
**
** This program does not really do much...it simply iterates until
** convergence or until iterations are exhausted.  It would not be
** necessary if it were possible to rewrite the PSI driver to be more
** general and allow non-crashing exits out of loops; somebody who
** knows how to do this should do it.
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
**
**
** Modification History:
**
** - Modified 10 February 1999 by C. David Sherrill -
** Added the ability to parse the orbital optimization log file (file14)
** so this information can be used to allow looser convergence on the 
** CI during early iterations.
**
*/

#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "setup_io.h"

#define MAX_COMMENT 10

void title(void);
void quote(void);
double calc_ci_conv(double);


FILE *infile, *outfile;
char *psi_file_prefix;


/* MAIN ROUTINE */

main(int argc, char *argv[])
{
  int converged = 0;
  int i, errcod = 0;
  int ncasiter = 0;            /* max cas iterations */
  char detci_string[80];       /* string containing system call for DETCI  */
  double ci_conv;              /* desired CI convergence 
                                  (changes dynamically during CAS opt)     */
  double scale_conv;           /* CI convergence threshold = 
                                  orbital gradient * scale_conv            */


  init_io(argc,argv);          /* open input and output files              */
  title();                     /* print program identification             */

  ncasiter = 30;
  errcod = ip_data("NCASITER","%d",&ncasiter,0);
  scale_conv = 0.01;
  errcod = ip_data("SCALE_CONV","%lf",&scale_conv,0);

  for (i=0; i<ncasiter && !converged; i++) {
    ci_conv = calc_ci_conv(scale_conv);

    if (ci_conv > 1.0E-7) {
      sprintf(detci_string, "detci --quiet -c %12.9lf\n", ci_conv);
    }
    else 
      sprintf(detci_string, "detci --quiet\n");

    check(!system("transqt --quiet"), "TRANSQT failed");
    check(!system(detci_string), "DETCI failed");
    check(!system("clag --quiet"), "CLAG failed");
    converged = system("detcas --quiet");
  }

  fprintf(outfile,"\n");
  fprintf(outfile,"*******************************************************\n");

  if (converged) 
    fprintf(outfile,"                  ORBITALS CONVERGED\n");
  else
    fprintf(outfile,"               ORBITALS DID NOT CONVERGE\n");

  if (converged) {
    system("rm -f detci_cfile.dat detci_sfile.dat");
    system("rm -f diis.dat orbs.dat thetas.dat");
  }

  quote();
  close_io();
  return(!converged);
}


/*
** title(): Function prints a program identification
*/
void title(void)
{
  fprintf(outfile,"\n");
  fprintf(outfile,"*******************************************************\n");
  fprintf(outfile,"                   D E T C A S M A N\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"                   C. David Sherrill\n") ;
  fprintf(outfile,"                    October 7 1998\n") ;
  fprintf(outfile,"*******************************************************\n");
  fprintf(outfile,"\n\n\n");
  fflush(outfile);
}


void quote(void)
{
  fprintf(outfile,"\n");
  fprintf(outfile,"                DETCAS MANAGER EXITING\n");
  fprintf(outfile,"*******************************************************\n");
  fprintf(outfile,"\n\n\n");
  fflush(outfile);
}


/*
** Read the current orbital convergence from file14
*/
double calc_ci_conv(double scale_conv)
{
  FILE *sumfile;
  char sumfile_name[] = "file14.dat";
  char comment[MAX_COMMENT];
  int i, entries, iter, nind;
  double scaled_rmsgrad, rmsgrad, energy_last;
  double tval;

  sumfile = fopen(sumfile_name, "r");

  if (sumfile == NULL) {
    return(1.0E-3);
  }

  if (fscanf(sumfile, "%d", &entries) != 1) {
    fprintf(outfile,"detcasman: Trouble reading num entries in file %s\n",
            sumfile_name);
    fclose(sumfile);
    return(1.0E-3);
  }

  for (i=0; i<entries; i++) {
    fscanf(sumfile, "%d %d %lf %lf %lf %s", &iter, &nind, &scaled_rmsgrad,
           &rmsgrad, &energy_last, comment);
  }
  fclose(sumfile);

  tval = (scaled_rmsgrad < rmsgrad) ? scaled_rmsgrad : rmsgrad;
  tval *= scale_conv;

  return(tval);
}

