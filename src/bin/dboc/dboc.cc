/*
**
**  DBOC: Driver program for computing the Diagonal Born-Oppenheimer Correction
**
*/

#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern "C" {
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/slaterdset.h>
#include <psifiles.h>
}
#include <libbasis/basisset.h>
#include <libbasis/overlap.h>
#include <masses.h>
#include <physconst.h>
#include "defines.h"
#include "molecule.h"
#include "moinfo.h"
#include "params.h"
#include "hfwfn.h"

#if defined HAVE_DECL_SETENV && !HAVE_DECL_SETENV
  extern int setenv(const char *, const char *, int);
#endif

/* Function prototypes */
static void init_io(int argc, char *argv[]);
static void exit_io();
void done(const char * message);
static void parsing();
static void read_molecule();
static double eval_dboc();
extern void print_intro();
extern void print_params();
extern "C" char *gprgid();
extern void setup_geoms();
extern double eval_derwfn_overlap();
extern void read_moinfo();

/*--- Global structures ---*/
FILE *infile, *outfile;
char *psi_file_prefix;
Molecule_t Molecule;
MOInfo_t MOInfo;
Params_t Params;

BasisSet* BasisSets[MAX_NUM_DISP];  // Array of pointers to basis set objects with current coordinate displaced by +delta, -delta, +2delta, and -2delta, respectively
HFWavefunction* HFVectors[MAX_NUM_DISP];  // Array of pointers to HF wavefunctions for displacements by +delta, -delta, +2delta, and -2delta, respectively
char* CI_Vector_Labels[MAX_NUM_DISP] = {
  "R0-delta CI Vector",
  "R0+delta CI Vector",
  "R0-2delta CI Vector",
  "R0+2delta CI Vector"
  };

const int MAX_GEOM_STRING=20;

int main(int argc, char *argv[])
{
  int i,j;
  int natom, num, junk;
  double **geom, *zvals;
  char *geom_string;
  FILE *geometry;

  init_io(argc,argv);
  parsing();
  print_intro();
  read_molecule();
  print_params();
  read_moinfo();
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

  Params.delta = 0.0005;
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
      if (Params.coords[coord].index < 0 || Params.coords[coord].index >= Molecule.natom*3)
	done("Keyword COORDS contains an out-of-bounds index");
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
  errcod = ip_data(":DBOC:DISP_PER_COORD","%d",&Params.disp_per_coord,0);
  if (Params.disp_per_coord != 2 && Params.disp_per_coord != 4)
    throw std::runtime_error("dboc.cc:parsing() -- disp_per_coord must either be 2 or 4");
}

/*--- Open chkpt file and grab molecule info ---*/
void read_molecule()
{
  chkpt_init(PSIO_OPEN_OLD);
  Molecule.natom = chkpt_rd_natom();
  Molecule.geom = chkpt_rd_geom();
  Molecule.zvals = chkpt_rd_zvals();
  chkpt_close();

  fprintf(outfile, "  -Reference Geometry:\n");
  for(int i=0; i < Molecule.natom; i++) {
    fprintf(outfile, "\n   %1.0f ", Molecule.zvals[i]);
    for(int j=0; j < 3; j++)
      fprintf(outfile, "%20.10f  ", Molecule.geom[i][j]);
  }
  fprintf(outfile, "\n\n");
}


// Returns an array of atomic masses, in a.u. (not in a.m.u.)
double* get_atomic_masses()
{
  double* atomic_mass;

  //
  // Convert isotope labels into atomic masses (in a.m.u.)
  //

  // Check number of atoms vs. number of isotope labels given in input
  if (Params.isotopes) {
    if (Molecule.natom != Params.nisotope)
      done("Number of atoms in molecule does not match the number of isotope labels in input file");
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

  return atomic_mass;
}


//
// This function will run psi for the first displaced point along the coordinate,
// save the basis set information, and save CI vector, if necessary
//
// NOTE: It's assumed that the first displacement is always by -delta
//
void run_psi_firstdisp(int disp)
{
  char *inputcmd = new char[80];
  sprintf(inputcmd,"input --geomdat %d",disp);
  int errcod = system(inputcmd);
  if (errcod) {
    done("input failed");
  }
  errcod = system("psi3 --dboc --noinput --messy");
  if (errcod) {
    done("psi3 failed");
  }
  delete[] inputcmd;

  //
  // Read in the "-delta" displaced basis
  //
  chkpt_init(PSIO_OPEN_OLD);
  BasisSets[MinusDelta] = new BasisSet(PSIF_CHKPT);
  // rotate the geometry back to the original frame and change shell centers
  double **rref_m = chkpt_rd_rref();
  double **geom_m = chkpt_rd_geom();
  chkpt_close();
  double **geom_m_ref = block_matrix(Molecule.natom,3);
  mmult(geom_m,0,rref_m,0,geom_m_ref,0,Molecule.natom,3,3,0);
  free_block(geom_m);
  free_block(rref_m);
  for(int a=0; a<Molecule.natom; a++)
    BasisSets[MinusDelta]->set_center(a,geom_m_ref[a]);
    
  // Read in the "-delta" displaced HF wavefunction
  HFVectors[MinusDelta] = new HFWavefunction();

  // For CI method rename the saved wave function
  if (!strcmp(Params.wfn,"DETCI") || !strcmp(Params.wfn,"DETCAS")) {
    SlaterDetVector *vec;
    slaterdetvector_read(PSIF_CIVECT,"CI vector",&vec);
    slaterdetvector_write(PSIF_CIVECT,CI_Vector_Labels[0],vec);
    slaterdetvector_delete_full(vec);
  }
}

void run_psi_otherdisp(int disp)
{
  char *inputcmd = new char[80];
  sprintf(inputcmd,"input --savemos --geomdat %d",disp);
  int errcod = system(inputcmd);
  if (errcod) {
    done("input failed");
  }
  errcod = system("psi3 --dboc --noinput --messy");
  if (errcod) {
    done("psi3 failed");
  }
  delete[] inputcmd;

  int disp_coord = (disp-1)%Params.disp_per_coord;

  // Read in the displaced HF wavefunction
  HFVectors[disp_coord] = new HFWavefunction();

  // For CI method rename the saved wave function
  if (!strcmp(Params.wfn,"DETCI") || !strcmp(Params.wfn,"DETCAS")) {
    SlaterDetVector *vec;
    slaterdetvector_read(PSIF_CIVECT,"CI vector",&vec);
    slaterdetvector_write(PSIF_CIVECT,CI_Vector_Labels[disp_coord],vec);
    slaterdetvector_delete_full(vec);
  }
}

//
// This function creates basis set objects for each displacement
//
void init_basissets(Params_t::Coord_t* coord)
{
  // BasisSetP1
  int atom = coord->index/3;
  int xyz = coord->index%3;
  double AplusD[3];
  for(int i=0; i<3; i++)
    AplusD[i] = BasisSets[MinusDelta]->get_center(atom, i);
  AplusD[xyz] += 2.0*Params.delta;
  BasisSets[PlusDelta] = new BasisSet(*BasisSets[MinusDelta]);
  BasisSets[PlusDelta]->set_center(atom,AplusD);
  
  if (Params.disp_per_coord == 4) {
    // BasisSetP2
    double Aplus2D[3];
    for(int i=0; i<3; i++)
      Aplus2D[i] = AplusD[i];
    Aplus2D[xyz] += Params.delta;
    BasisSets[Plus2Delta] = new BasisSet(*BasisSets[MinusDelta]);
    BasisSets[Plus2Delta]->set_center(atom,Aplus2D);

    // BasisSetM2
    double Aminus2D[3];
    for(int i=0; i<3; i++)
      Aminus2D[i] = Aplus2D[i];
    Aminus2D[xyz] -= 4.0*Params.delta;
    BasisSets[Minus2Delta] = new BasisSet(*BasisSets[MinusDelta]);
    BasisSets[Minus2Delta]->set_center(atom,Aminus2D);
  }
}

double eval_dboc()
{

  double* atomic_mass = get_atomic_masses();  
  const int ndisp = Params.disp_per_coord * Params.ncoord;
  double E_dboc = 0.0;

  for(int disp=1; disp<=ndisp;) {
    int current_coord = (disp-1)/Params.disp_per_coord;
    Params_t::Coord_t* coord = &(Params.coords[current_coord]);
    int atom = coord->index/3;
    int xyz = coord->index%3;

    run_psi_firstdisp(disp);
    disp++;
    for(int i=1;i<Params.disp_per_coord; i++) {
      run_psi_otherdisp(disp);
      disp++;
    }
    
    init_basissets(coord);

    // The - sign comes from the integration by parts
    double del2 = (-1.0)*eval_derwfn_overlap();
    int Z = (int)Molecule.zvals[atom];
    // mass of nucleus = atomic mass - mass of electrons
    double nuclear_mass = atomic_mass[atom]  - Z;
    double E_i = -del2/(2.0*nuclear_mass);
    // multiply by the degeneracy factor
    E_i *= coord->coeff;
    if (Params.print_lvl >= PrintLevels::print_params) {
      fprintf(outfile,"  <d2/dx2>                 = %25.15lf\n",del2);
      fprintf(outfile,"  DBOC contribution from cartesian degree of freedom %d = %20.10lf a.u.\n\n",
	      coord->index+1,E_i);
      fflush(outfile);
    }
    E_dboc += E_i;

    // For CI method purge the file with saved wave functions as the point group may change
    if (!strcmp(Params.wfn,"DETCI") || !strcmp(Params.wfn,"DETCAS")) {
      psio_open(PSIF_CIVECT,PSIO_OPEN_NEW);
      psio_close(PSIF_CIVECT,0);
    }
    // it's safe to do psiclean now
    int errcod = system("psiclean");
    if (errcod) {
      done("psiclean");
    }

    delete_basissets(coord);
    delete_hfwfns(coord);

  }

  return E_dboc;
}

static char *orig_psi_output_env;

void init_io(int argc, char *argv[])
{
  int i;
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  int errcod = psi_start(argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    throw std::runtime_error("init_io -- psi_start failed");
  ip_cwk_add(progid);
  tstart(outfile);
  psio_init();

  // Psi modules called by dboc should read from the same input file
  // set the value of PSI_INPUT for the duration of this run
  char* ifname = psi_ifname();
#if HAVE_PUTENV
  char* tmpstr1 = (char *) malloc(11+strlen(ifname));
  sprintf(tmpstr1, "PSI_INPUT=%s", ifname);
  putenv(tmpstr1);
#elif HAVE_SETENV
  setenv("PSI_INPUT",ifname,1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif

  // Psi modules called by dboc should write to a different output file
  // reset the value of PSI_OUTPUT for the duration of this run
  orig_psi_output_env = getenv("PSI_OUTPUT");
#if HAVE_PUTENV
  char* tmpstr2 = strdup("PSI_OUTPUT=dboc.findif.out");
  putenv(tmpstr2);
#elif HAVE_SETENV
  setenv("PSI_OUTPUT","dboc.findif.out",1);
#else
#error "Have neither putenv nor setenv. Something must be very broken on this system."
#endif

  free(progid);
}

void exit_io()
{
  psio_done();
  tstop(outfile);
  psi_stop();
}

void done(const char *message)
{
  char* errmsg;
  errmsg = (char *) new char[strlen(message)+7];
  sprintf(errmsg,"DBOC: %s",message);
  exit_io();
  throw std::runtime_error(errmsg);
  delete[] errmsg;
}


extern "C" char *gprgid()
{
   char *prgid = "DBOC";

   return(prgid);
}
