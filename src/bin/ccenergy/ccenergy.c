/*
**  CCENERGY: Program to calculate coupled cluster energies.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <dpd.h>
#include <file30.h>
#include "globals.h"

/* Function prototypes */
void init_io(void);
void title(void);
void get_moinfo(void);
void get_params(void);
void init_amps(void);
void tau_build(void);
void taut_build(void);
double energy(void);
void sort_amps(void);
void Fae_build(void);
void Fmi_build(void);
void Fme_build(void);
void t1_build(void);
void Wmnij_build(void);
void Z_build(void);
void Y_build(void);
void X_build(void);
void Wmbej_build(void);
void t2_build(void);
void tsave(void);
int converged(void);
double diagnostic(void);
double d1diag(void);
void exit_io(void);
void cleanup(void);
void update(void);
void diis(int iter);
void ccdump(void);

int main(int argc, char *argv[])
{
  int done=0;
  int i, natom;
  double **geom, *zvals;
  FILE *efile;

  moinfo.iter=0;
  
  init_io();
  title();
  get_moinfo();
  get_params();
  dpd_init(moinfo.nirreps, params.memory, 2, moinfo.occpi, moinfo.occ_sym,
	   moinfo.virtpi, moinfo.vir_sym);
  init_amps();
  tau_build();
  taut_build();
  fprintf(outfile, "\t                     Solving CCSD Equations\n");
  fprintf(outfile, "\t                     ----------------------\n");
  fprintf(outfile, "\tIter             Energy               RMS       T1Diag      D1Diag\n");
  fprintf(outfile, "\t----     ---------------------     --------   ----------  ----------\n");
  moinfo.ecc = energy();
  moinfo.t1diag = diagnostic();
  moinfo.d1diag = d1diag();
  update();
  for(moinfo.iter=1; moinfo.iter <= params.maxiter; moinfo.iter++) {
      sort_amps();
      Y_build();
      X_build();
      Wmbej_build(); 
      Fme_build(); Fae_build(); Fmi_build();
      t1_build();
      Wmnij_build(); Z_build();
      t2_build();
      if(converged()) {
	  done = 1;  /* Boolean for convergence */
	  tsave();
	  tau_build(); taut_build();
	  moinfo.ecc = energy();
	  moinfo.t1diag = diagnostic();
	  moinfo.d1diag = d1diag();
	  sort_amps();  /* For upcoming calculations */
	  update();
	  fprintf(outfile, "\n\tIterations converged.\n");
          fprintf(outfile, "\n");
          d1diag_print();
	  fflush(outfile);
	  break;
	}
      diis(moinfo.iter); 
      tsave();
      tau_build(); taut_build();
      moinfo.ecc = energy();
      moinfo.t1diag = diagnostic();
      moinfo.d1diag = d1diag();
      update();
    }
  fprintf(outfile, "\n");
  if(!done) {
     fprintf(outfile, "\t ** Wave function not converged to %2.1e ** \n",
	     params.convergence);
     fflush(outfile);
     dpd_close();
     cleanup();
     exit_io();
     exit(1);
    }

  fprintf(outfile, "\tSCF energy       (file30)  = %20.15f\n", moinfo.escf);
  fprintf(outfile, "\tReference energy (file100) = %20.15f\n", moinfo.eref);
  fprintf(outfile, "\tCCSD correlation energy    = %20.15f\n", moinfo.ecc);
  fprintf(outfile, "\tTotal CCSD energy          = %20.15f\n", 
          moinfo.eref + moinfo.ecc);
  fprintf(outfile, "\n");

/*  Curt had these print statements in the code.  They are not
**  necessary now, but I just comment them out for the time being. - MLL
**  12-3-99
    {
    struct dpdbuf tIJAB, tijab, tIjAb;
    dpd_buf_init(&tIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
    dpd_buf_print(&tIJAB, outfile);
    dpd_buf_close(&tIJAB);
    dpd_buf_init(&tijab, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
    dpd_buf_print(&tijab, outfile);
    dpd_buf_close(&tijab);
    dpd_buf_init(&tIjAb, CC_TAMPS, 2, 7, 2, 7, 0, "tIjAb", 0, outfile);
    dpd_buf_print(&tIjAb, outfile);
    dpd_buf_close(&tIjAb);
  }
*/
/* Stop here */



  /* Write pertinent data to energy.dat for Dr. Yamaguchi */
  file30_init();
  natom = file30_rd_natom();
  geom = file30_rd_geom();
  zvals = file30_rd_zvals();
  file30_close();
  ffile(&efile, "energy.dat",1);
  fprintf(efile, "*\n");
  for(i=0; i < natom; i++) 
      fprintf(efile, " %4d   %5.2f     %13.10f    %13.10f    %13.10f\n",
	      i+1, zvals[i], geom[i][0], geom[i][1], geom[i][2]);
  free_matrix(geom, natom);  free(zvals);
  fprintf(efile, "SCF(30)   %22.12f\n", moinfo.escf);
  fprintf(efile, "REF(100)  %22.12f\n", moinfo.eref);
  fprintf(efile, "CCSD      %22.12f\n", (moinfo.ecc+moinfo.eref));
  fclose(efile);
  
  dpd_close();
  cleanup(); 
  exit_io();
  exit(0);
}

void init_io(void)
{
  int i;
  char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(progid);

  free(progid);

  psio_init();
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCENERGY        *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);
  psio_done();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "CCENERGY";

   return(prgid);
}
