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
#include <libfile30/file30.h>
#include <libpsio/psio.h>
}
#include "molecule.h"
#include "params.h"

/* Function prototypes */
static void init_io();
static void exit_io();
void done(const char * message);
static void parsing();
static void read_molecule();
extern "C" char *gprgid();
extern void setup_geoms();
extern void eval_dboc_rhf();

/*--- Global structures ---*/
FILE *infile, *outfile;
Molecule_t Molecule;
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
  read_molecule();
  setup_geoms();
  if (!strcmp(Params.wfn,"SCF")) {
    if (Params.reftype == Params_t::rhf) {
      eval_dboc_rhf();
    }
    else
      done("Only RHF reference is supported at the moment");
  }
  exit_io();
  exit(0);
}


/*--- parsing ---*/
void parsing()
{
  int errcod;
  
  errcod = ip_string("WFN",&Params.wfn,0);
  if (errcod != IPE_OK)
    done("Keyword WFN is not found");

  char *reftype;
  errcod = ip_string("REFERENCE",&reftype,0);
  if (errcod != IPE_OK)
    done("Keyword REFERENCE is not found");
  else if (!strcmp(reftype,"RHF"))
    Params.reftype = Params_t::rhf;
  else
    done("Only RHF reference is supported at the moment");
  delete[] reftype;

  Params.delta = 0.00001;
  Params.disp_per_coord = 2;
}

/*--- Open file30 and grab molecule info ---*/
void read_molecule()
{
  file30_init();
  Molecule.natom = file30_rd_natom();
  Molecule.geom = file30_rd_geom();
  Molecule.zvals = file30_rd_zvals();
  int nirreps = file30_rd_nirreps();
  if (nirreps != 1)
    done("DBOC computations currently possible only in C1 symmetry");
  file30_close();

  fprintf(outfile, "Reference Geometry:\n");
  for(int i=0; i < Molecule.natom; i++) {
    fprintf(outfile, "\n %1.0f ", Molecule.zvals[i]);
    for(int j=0; j < 3; j++)
      fprintf(outfile, "%20.10f  ", Molecule.geom[i][j]);
  }
  fprintf(outfile, "\n\n");
}

void init_io(void)
{
  int i;
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
