/*
**
**  DBOC: Driver program for computing the Diagonal Born-Oppenheimer Correction
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern "C" {
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/slaterd.h>
#include <psifiles.h>
}
#include <masses.h>
#include <physconst.h>
#include "molecule.h"
#include "moinfo.h"
#include "params.h"

/* Function prototypes */
static void init_io();
static void exit_io();
void done(const char * message);
static void parsing();
static void read_chkpt();
static double eval_dboc();
extern void print_intro();
extern void print_params();
extern "C" char *gprgid();
extern void setup_geoms();
extern double eval_derwfn_overlap();

/*--- Global structures ---*/
FILE *infile, *outfile;
Molecule_t Molecule;
MOInfo_t MOInfo;
Params_t Params;

#define MAX_GEOM_STRING 20

int main(int argc, char *argv[])
{
  int i,j;
  int natom, num, junk;
  double **geom, *zvals;
  char *geom_string;
  FILE *geometry;

  init_io();
  parsing();
  print_intro();
  print_params();
  read_chkpt();
  setup_geoms();
  double E_dboc = eval_dboc();
  fprintf(outfile,"  E(DBOC) = %25.15lf a.u.\n",E_dboc);
  fprintf(outfile,"  E(DBOC) = %25.5lf cm^{-1}\n\n",E_dboc*_hartree2wavenumbers);
  exit_io();
  exit(0);
}


/*--- parsing ---*/
void parsing()
{
  int errcod;
  
  errcod = ip_string("LABEL",&Params.label,0);
  if (errcod != IPE_OK) {
    Params.label = new char[1];
    Params.label[0] = '\0';
  }

  errcod = ip_string("WFN",&Params.wfn,0);
  if (errcod != IPE_OK)
    done("Keyword WFN is not found");

  char *reftype;
  errcod = ip_string("REFERENCE",&reftype,0);
  if (errcod != IPE_OK)
    done("Keyword REFERENCE is not found");
  else if (!strcmp(reftype,"RHF"))
    Params.reftype = Params_t::rhf;
  else if (!strcmp(reftype,"ROHF"))
    Params.reftype = Params_t::rohf;
  else if (!strcmp(reftype,"UHF"))
    Params.reftype = Params_t::uhf;
  else
    done("This HF reference is not supported at the moment");
  delete[] reftype;

  Params.delta = 0.0001;
  errcod = ip_data(":DBOC:DISPLACEMENT","%lf",&Params.delta,0);

  Params.print_lvl = 1;
  errcod = ip_data("PRINT","%d",&Params.print_lvl,0);

  Params.disp_per_coord = 2;
}

/*--- Open file30 and grab molecule info ---*/
void read_chkpt()
{
  chkpt_init();

#if 0
  Molecule.natom = file30_rd_natom();
  Molecule.geom = file30_rd_geom();
  Molecule.zvals = file30_rd_zvals();
  int nirreps = file30_rd_nirreps();
  if (nirreps != 1)
    done("DBOC computations currently possible only in C1 symmetry");

  MOInfo.num_so = file30_rd_nso();
  MOInfo.num_mo = file30_rd_nmo();
  int *clsdpi = file30_rd_clsdpi();
  MOInfo.ndocc = clsdpi[0];
  delete[] clsdpi;
  int *openpi = file30_rd_openpi();
  MOInfo.nsocc = openpi[0];
  delete[] openpi;
  MOInfo.nalpha = MOInfo.ndocc + MOInfo.nsocc;
  MOInfo.nbeta = MOInfo.ndocc;
#else
  Molecule.natom = chkpt_rd_natom();
  Molecule.geom = chkpt_rd_geom();
  Molecule.zvals = chkpt_rd_zvals();
  int nirreps = chkpt_rd_nirreps();
  if (nirreps != 1)
    done("DBOC computations currently possible only in C1 symmetry");

  MOInfo.num_so = chkpt_rd_nso();
  MOInfo.num_mo = chkpt_rd_nmo();
  int *clsdpi = chkpt_rd_clsdpi();
  MOInfo.ndocc = clsdpi[0];
  delete[] clsdpi;
  int *openpi = chkpt_rd_openpi();
  MOInfo.nsocc = openpi[0];
  delete[] openpi;
#endif

  chkpt_close();

  fprintf(outfile, "  -Reference Geometry:\n");
  for(int i=0; i < Molecule.natom; i++) {
    fprintf(outfile, "\n   %1.0f ", Molecule.zvals[i]);
    for(int j=0; j < 3; j++)
      fprintf(outfile, "%20.10f  ", Molecule.geom[i][j]);
  }
  fprintf(outfile, "\n\n");
}

double eval_dboc()
{
  const int ndisp = Params.disp_per_coord * Molecule.natom * 3;
  double E_dboc = 0.0;

  for(int disp=1; disp<=ndisp; disp++) {
    char *inputcmd = new char[80];
    int atom = (disp-1)/6;
    sprintf(inputcmd,"input --geomdat %d --noreorient --nocomshift",disp);
    int errcod = system(inputcmd);
    if (errcod) {
      done("input failed");
    }
    errcod = system("psi");
    if (errcod) {
      done("psi failed");
    }
    disp++;

    // For CI method rename the saved wave function
    if (!strcmp(Params.wfn,"DETCI")) {
      SlaterDetVector *vec;
      slaterdetvector_read(PSIF_CIVECT,"CI vector",&vec);
      slaterdetvector_write(PSIF_CIVECT,"Old CI vector",vec);
      slaterdetvector_delete_full(vec);
    }

    sprintf(inputcmd,"input --savemos --geomdat %d --noreorient --nocomshift",disp);
    errcod = system(inputcmd);
    if (errcod) {
      done("input failed");
    }
    errcod = system("psi");
    if (errcod) {
      done("psi failed");
    }
    delete[] inputcmd;

    double S = eval_derwfn_overlap();
    double del2 = (1.0-S)/(2.0*Params.delta*Params.delta);
    double E_i = del2*_au2amu/(2.0*an2masses[(int)Molecule.zvals[atom]]);
    if (Params.print_lvl > PrintLevels::print_intro) {
      fprintf(outfile,"  +- wave function overlap = %25.15lf\n",S);
      fprintf(outfile,"  DBOC contribution from cartesian degree of freedom %d = %20.10lf a.u.\n\n",
	      (disp-1)/2+1,E_i);
      fflush(outfile);
    }
    E_dboc += E_i;
  }

  return E_dboc;
}

void init_io(void)
{
  int i;
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"dboc.out",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(progid);
  psio_init();

  free(progid);
}

void exit_io()
{
  int i;

  psio_done();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

void done(const char *message)
{
  fprintf(outfile,"%s\n",message);
  fprintf(stderr,"DBOC: %s\n",message);
  exit_io();
  abort();
}


extern "C" char *gprgid()
{
   char *prgid = "DBOC";

   return(prgid);
}
