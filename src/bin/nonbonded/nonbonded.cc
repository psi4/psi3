/*
** NONBONDED
**
** Evaluate nonbonded interactions such as van der Waals terms and
** point-charge electrostatics
**
** C. David Sherrill
** January 2008
*/


/*! 
** \file
** \ingroup NONBONDED
** \brief Evaluate empirical non-bonded interactions
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cctype>               // for toupper()
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <cstring>
#include <physconst.h>
#include <psifiles.h>
#include <masses.h>

#include "globals.h"
#include "nonbonded.h"

namespace psi { namespace nonbonded {

/* function prototypes this module */
void start_io(int argc, char *argv[]);
void stop_io(void);
void print_intro(void);
void compute_R(int natom, double **geom, double *R);
double compute_ddisp(int natom, double *R, double *AN, double s6, double d,
  int damp_flag);
int label2an(char *label);

}}; // close namespace decl

using namespace psi::nonbonded;

int main(int argc, char *argv[]) {

  double **geom;                  /* geometry matrix (cols are x,y,z)   */
  double *AC;                     /* atomic charges                     */
  double *AN;                     /* atomic numbers                     */
  double *R;                      /* distance matrix (lwr triangle)     */
  int natom;                      /* number of atoms                    */
  int have_partial_charges=0;     /* flag for partial charges available */
  double s6=1.0;                  /* global scaling parameter (Grimme),
                                     0.75 PBE, 1.2 BLYP, 1.05 B-P86,
                                     1.0 TPSS, 1.05 B3LYP               */
  double d=20.0;                  /* exponent for damping term (Grimme) */
  double energy_dd;               /* damped dispersion energy (J/mol)   */
  double energy_dd_hartree;       /* in hartree                         */
  double e_scf;                   /* Hartree-Fock energy (hartree)      */
  int damp_flag;                  /* damp the dispersion? 1=yes,0=no    */
  int num_array;                  /* count elements in parsed array     */
  int errcod;                     /* input parsing error code           */
  double tval;                    /* temp var                           */
  char tmpstr[100];               /* temp string for parsing input      */
  int i,j;


  start_io(argc,argv);
  print_intro();
  geom = chkpt_rd_geom();
  natom = chkpt_rd_natom();
  AN = chkpt_rd_zvals();
  if (ip_exist("PARTIAL_CHARGES",0)) {
    AC = init_array(natom);
    errcod = ip_double_array("PARTIAL_CHARGES",AC,natom);
    have_partial_charges = 1;
  }
  errcod = ip_data("S6","%lf",&s6,0);
  errcod = ip_data("D_DMP","%lf",&d,0);
  damp_flag = 1;
  errcod = ip_boolean("DAMP",&(damp_flag),0);

  /* parse overridden vdW radii */
  if (ip_exist("VDW_RADII",0)) {
    errcod = ip_count("VDW_RADII", &num_array, 0);
    if (errcod != IPE_OK) {
      fprintf(stderr,"ERROR: %s\n", ip_error_message(errcod));
    }  
    else {
      for (i=0; i<num_array; i+=2) {
        errcod = ip_data("VDW_RADII", "%s", tmpstr, 1, i); 
        errcod = ip_data("VDW_RADII", "%lf", &tval, 1, i+1); 
        if ((j = label2an(tmpstr)) != -1) {
          vdw_radii_grimme[j] = tval;
          fprintf(outfile, 
            "  Using custom value of %6.3lf for vdW radius of element %s\n",
            tval, tmpstr);
        }
      }
    }
  } /* end parsing overridden vdW radii */

  /* parse overridden C6 coefficients */
  if (ip_exist("C6",0)) {
    errcod = ip_count("C6", &num_array, 0);
    if (errcod != IPE_OK) {
      fprintf(stderr,"ERROR: %s\n", ip_error_message(errcod));
    }  
    else {
      for (i=0; i<num_array; i+=2) {
        errcod = ip_data("C6", "%s", tmpstr, 1, i); 
        errcod = ip_data("C6", "%lf", &tval, 1, i+1); 
        if ((j = label2an(tmpstr)) != -1) {
          vdw_C6_grimme[j] = tval;
          fprintf(outfile, 
            "  Using custom value of %6.3lf for C6 for element %s\n",
            tval, tmpstr);
        }
      }
    }
  } /* end parsing overridden C6 coefficients */
 
  fprintf(outfile, "\n");
  fprintf(outfile, "  Universal scaling coefficient s6 = %6.4lf\n", s6);
  fprintf(outfile, "  Universal damping exponent d     = %6.4lf\n", d);
  fprintf(outfile, "  Damp dispersion                  = %6s\n", 
    (damp_flag==1) ? "yes" : "no");

  R = init_array((natom*(natom+1))/2);
  compute_R(natom, geom, R);
  energy_dd = compute_ddisp(natom, R, AN, s6, d, damp_flag);
  energy_dd_hartree = energy_dd / (_na * _hartree2J);

  fprintf(outfile, "\n");
  fprintf(outfile, "  Damped dispersion energy  = %14.9lf hartree ", 
    energy_dd_hartree);
  fprintf(outfile, "(%14.9lf kcal/mol)\n", (energy_dd / 4184));
  fflush(outfile);

  /* get the SCF energy so we can print HF+D */
  if (psio_tocscan(PSIF_CHKPT, "SCF energy") != NULL) {
    e_scf = chkpt_rd_escf();
    fprintf(outfile, "  Hartree-Fock energy       = %14.9lf hartree\n", e_scf);
    fprintf(outfile, "  Hartree-Fock + dispersion = %14.9lf hartree\n", 
      e_scf + energy_dd_hartree);
  }

  fprintf(outfile, "\n");

  /* clean up */
  free(AN);
  free(R);
  if (have_partial_charges) free(AC);
  /* free geom --- how exactly is it allocated now? */
  stop_io();
}

extern "C" { char *gprgid(void) { char *prgid = "TRANSQT"; return (prgid); } }

namespace psi { namespace nonbonded {

/*!
** start_io(): Initiate PSI IO routines
**
** \param argc = number of command-line arguments
** \param argv = the command-line arguments
**
** \ingroup NONBONDED
*/
void start_io(int argc, char *argv[])
{
  int errcod;
  
  errcod = psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    abort();
  ip_cwk_add(":NONBONDED");
  tstart(outfile);
  psio_init(); psio_ipv1_config();
  chkpt_init(PSIO_OPEN_OLD);
   
  return;
}


/*!
** stop_io(): Shut down PSI IO
** 
** Returns: none
**
** \ingroup NONBONDED
*/
void stop_io(void)
{
  tstop(outfile);
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);
}
 


/*!
** print_intro(): Print an intro for the module
**
** Returns: none
**
** \ingroup NONBONDED
*/
void print_intro(void)
{
 fprintf(outfile,"             ---------------------------------------\n");
 fprintf(outfile,"                            NONBONDED                  \n");
 fprintf(outfile,"               Evaluate empirical non-bonded terms     \n");
 fprintf(outfile,"                        C. David Sherrill              \n");
 fprintf(outfile,"             ---------------------------------------\n");
 fprintf(outfile,"\n");
}


/*!
** compute_R(): Compute the distances between each pair of atoms.  
** Convert to Angstrom because those are the units Grimme's parameters 
** are in.
**
** \param natom = number of atoms
** \param geom  = geometry matrix (cols are x, y, z; rows are atoms) 
** \param R     = matrix of distances (packed lower triangle, ij indexing)
**
** Returns: none
**
** \ingroup NONBONDED
*/
void compute_R(int natom, double **geom, double *R)
{
  int i, j, ij;
  double dx, dy, dz, dist;

  for (i=0,ij=0; i<natom; i++) {
    for (j=0; j<i; j++,ij++) {
      dx = geom[i][0] - geom[j][0];
      dy = geom[i][1] - geom[j][1];
      dz = geom[i][2] - geom[j][2];
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      R[ij] = dist * _bohr2angstroms;
    }
  }
 
}



/*!
** compute_ddisp(): Compute damped dispersion terms 
** (ala S. Grimme, J. Comput. Chem. 27, 1787, 2006)
**
** \param natom      = number of atoms
** \param R          = matrix of interatomic distances (packed lower triangle)
** \param AN         = atomic number for each atom
** \param s6         = s6 universal dispersion scaling parameter
** \param d          = d universal dispersion damping exponent
** \param damp_flag  = 1 to damp dispersion, 0 otherwise (should be 1
**                     except for testing)
**
** Returns: the dispersion energy, damped if damp_flag=1
**
** \ingroup NONBONDED
*/
double compute_ddisp(int natom, double *R, double *AN, double s6, double d,
  int damp_flag)
{
  int i, j, ij, Zi, Zj;
  double energy=0.0, tval, r, r_i, r_j, r_vdw;
  double C6i, C6j, fdmp; 

  /* loop over unique pairs of atoms and evaluate the C6 dispersion terms */
  for (i=0,ij=0; i<natom; i++) {
    Zi = (int) AN[i];
    C6i = vdw_C6_grimme[Zi];
    r_i = vdw_radii_grimme[Zi];
    for (j=0; j<i; j++,ij++) {
      Zj = (int) AN[j];
      r = R[ij];
      C6j = vdw_C6_grimme[Zj];
      r_j = vdw_radii_grimme[Zj];
      r_vdw = r_i + r_j;
      tval = sqrt(C6i * C6j) * 1000000.0; /* nm^6 -> Ang^6 */
      tval = tval/pow(r, 6.0);
      printf("Undamped Dispersion term: %14.9lf\n", tval / 4184);
      if (damp_flag) {
        fdmp = 1.0 / (1.0 + exp(-d * (r / r_vdw - 1.0))); 
        printf("Damping factor          : %14.9lf\n", fdmp);
        tval *= fdmp;
      }
      energy += tval;
    }
  } 

  energy = -s6 * energy; /* in J mol^-1 */
  return energy;
}


/*!
** label2an(): This function returns the atomic number corresponding to a 
** given mass label, or -1 if not found.
**
** \param label = Atom label we're trying to match (case insensitive)
**
** Returns: atomic number of matched label (or -1 if not found)
**
** Depends on atomic_labels in masses.h
**
** C. David Sherrill
** July 1999
**
** \ingroup NONBONDED
*/
int label2an(char *label)
{
  int i, j, k;
  char p, q;

  for (i=0; i<LAST_ATOMIC_INDEX; i++) {
    k = strlen(label);
    for (j=0; j < k; j++) {
       p = label[j];
       q = atomic_labels[i][j];
       if (toupper(p) != toupper(q)) break;
    }
    if (j == k) return(i);
  }

  /* couldn't find the label */
  return(-1);

} // end label2an()


}} // close off psi::nonbonded
