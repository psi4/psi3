/*
** READ30: A program to read data from PSI's file30.
*/

#include <stdio.h>
#include "file30.h"
#include <libciomr.h>
#include <ip_libv1.h>

FILE *infile, *outfile;
void init_io(void);
void title(void);
void exit_io(void);

void main(void)
{
  int i;
  int *opi, *cpi, *spi;
  double **geom,**scf, **usotao;
  double *zvals, *evals;
  char **irr_labs;
  
  init_io();
  title();

  file30_init();

  fprintf(outfile, "\tFile30 is %llu bytes long.\n", flen(30));

  fprintf(outfile,"\n\tLabel:       %s\n", file30_rd_label());
  fprintf(outfile,  "\tNum. atoms:  %d\n", file30_rd_natom());
  fprintf(outfile,  "\tNum. MOs:    %d\n", file30_rd_nmo());
  fprintf(outfile,  "\tNum. AOs:    %d\n", file30_rd_nao());
  fprintf(outfile,  "\tNum. irreps: %d\n", file30_rd_nirreps());
  fprintf(outfile,  "\tNum. irreps in HF: %d\n", file30_rd_nsymhf());
  fprintf(outfile,  "\tNum. open shells:  %d\n", file30_rd_iopen());
  fprintf(outfile,  "\tNum. allowed HF coeffs: %d\n", file30_rd_mxcoef());

  fprintf(outfile,"\n\tMolecular geometry:\n");
  fprintf(outfile,  "\t-------------------\n");
  geom = file30_rd_geom();
  zvals = file30_rd_zvals();

  for(i=0; i < file30_rd_natom(); i++) {
      fprintf(outfile, "     %1.0f  %15.10f %15.10f %15.10f\n", zvals[i],
	      geom[i][0], geom[i][1], geom[i][2]);
    }

  free_matrix(geom,file30_rd_natom());
  free(zvals);

  fprintf(outfile,"\n\tHartree-Fock energy:      %20.10f\n",file30_rd_escf());
  fprintf(outfile,  "\tNuclear repulsion energy: %20.10f\n",file30_rd_enuc());
  fprintf(outfile,  "\tCorrelation energy:       %20.10f\n",file30_rd_ecorr());

  opi = file30_rd_orbspi();
  cpi = file30_rd_clsdpi();
  spi = file30_rd_openpi();
  irr_labs = file30_rd_irr_labs();
  fprintf(outfile,"\n\tLabel\t# MOs\t# DOCC\t# SOCC\n");
  fprintf(outfile,  "\t-----\t-----\t------\t------\n");
  for(i=0; i < file30_rd_nirreps(); i++)
      fprintf(outfile,"\t %s\t   %d\t    %d\t    %d\n",
		  irr_labs[i],opi[i],cpi[i],spi[i]);

  fprintf(outfile,"\n\tMolecular orbital coefficients and energies:\n");
  fprintf(outfile,  "\t--------------------------------------------\n");
  scf = file30_rd_scf();
  evals = file30_rd_evals();
  eivout(scf,evals,file30_rd_nao(),file30_rd_nmo(),outfile);
  free_matrix(scf,file30_rd_nao());
  free(evals,file30_rd_nmo());

  file30_close();

  exit_io();
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*         READ30         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

char *gprgid()
{
   char *prgid = "READ30"; return(prgid);
}

void init_io(void)
{
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":READ30");
}

void exit_io(void)
{
  tstop(outfile);
  ip_done();
}
