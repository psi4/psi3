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

  int coords_given = ip_exist(":DBOC:COORDS",0);
  if (coords_given) {
    Params.ncoord = 0;
    ip_count(":DBOC:COORDS",&Params.ncoord,0);
    if (Params.ncoord == 0)
      done("Keyword COORDS should be a vector of 2-element vectors");
    Params.coords = new Params_t::Coord_t[Params.ncoord];
    for(int coord=0; coord<Params.ncoord; coord++) {
      int nelem = 0;
      ip_count(":DBOC:COORDS",&nelem,1,coord);
      if (nelem != 2)
	done("Keyword COORDS should be a vector of 2-element vectors");
      errcod = ip_data(":DBOC:COORDS","%d",&Params.coords[coord].index,2,coord,0);
      errcod = ip_data(":DBOC:COORDS","%lf",&Params.coords[coord].coeff,2,coord,1);
    }
  }
  else {
    Params.ncoord = 0;
    Params.coords = NULL;
  }

  int isotopes_given = ip_exist(":DBOC:ISOTOPES",0);
  if (isotopes_given) {
    Params.nisotope = 0;
    ip_count(":DBOC:ISOTOPES",&Params.nisotope,0);
    if (Params.nisotope == 0)
      done("Keyword ISOTOPES should be a vector of num_atoms elements");
    Params.isotopes = new char*[Params.nisotope];
    for(int atom=0; atom<Params.nisotope; atom++) {
      ip_string(":DBOC:ISOTOPES",&(Params.isotopes[atom]),1,atom);
    }
  }
  else {
    Params.nisotope = 0;
    Params.isotopes = NULL;
  }

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
  MOInfo.nalpha = MOInfo.ndocc + MOInfo.nsocc;
  MOInfo.nbeta = MOInfo.ndocc;
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
  double* atomic_mass;

  //
  // Convert isotope labels into atomic masses (in a.m.u.)
  //
  // Check number of atoms vs. number of isotope labels given in input
  if (Params.isotopes) {
    if (Molecule.natom != Params.nisotope)
      done("Number of atoms in molecule does not match the number of isotope labels in input.dat");
    else {
      atomic_mass = new double[Molecule.natom];
      for(int atom=0; atom<Molecule.natom; atom++) {
	char* isotope_label = Params.isotopes[atom];
	int label;
	for(label=0; label<=LAST_MASS_INDEX; label++) {
	  if (!strcmp(mass_labels[label],isotope_label))
	    break;
	}
	atomic_mass[atom] = atomic_masses[label];
      }
    }
  }
  else {
    atomic_mass = new double[Molecule.natom];
    for(int atom=0; atom<Molecule.natom; atom++)
      atomic_mass[atom] = an2masses[(int)Molecule.zvals[atom]];
  }
  
  //
  // Convert atomic masses to a.u.
  //
  for(int atom=0; atom<Molecule.natom; atom++)
    atomic_mass[atom] /= _au2amu;

  const int ndisp = Params.disp_per_coord * Params.ncoord;
  double E_dboc = 0.0;

  for(int disp=1; disp<=ndisp; disp++) {
    char *inputcmd = new char[80];
    Params_t::Coord_t* coord = &(Params.coords[(disp-1)/2]);
    int atom = coord->index/3;
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
    double del2 = (S-1.0)/(2.0*Params.delta*Params.delta);
    int Z = (int)Molecule.zvals[atom];
    // mass of nucleus = atomic mass - mass of electrons
    double nuclear_mass = atomic_mass[atom]  - Z;
    double E_i = -del2/(2.0*nuclear_mass);
    // multiply by the degeneracy factor
    E_i *= coord->coeff;
    if (Params.print_lvl > PrintLevels::print_intro) {
      fprintf(outfile,"  +- wave function overlap = %25.15lf\n",S);
      fprintf(outfile,"  <d2/dx2>                 = %25.15lf\n",del2);
      fprintf(outfile,"  DBOC contribution from cartesian degree of freedom %d = %20.10lf a.u.\n\n",
	      coord->index+1,E_i);
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
