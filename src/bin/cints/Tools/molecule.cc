/*! \file molecule.cc
    \ingroup CINTS
    \brief Enter brief description of file here
*/
#include<cstdio>
#include<cmath>
#include<cstdlib>
/* isnan and isinf are defined here in IBM C/C++ compilers */
#ifdef __IBMCPP__
#include<math.h>
#endif
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include<libchkpt/chkpt.hpp>
#include<libint/libint.h>
#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static void get_geometry(void);

void init_molecule()
{
  Molecule.label = chkpt_rd_label();
  Molecule.num_atoms = chkpt_rd_natom();
  Molecule.Rref = chkpt_rd_rref();
  get_geometry();

  return;
}


void cleanup_molecule()
{
  free(Molecule.centers);
  free(Molecule.Rref);

  return;
}

void get_geometry()
{
   int i;
   double *Z;  /* nuclear charges */
   double **g; /* cartesian geometry */

   Molecule.centers = (struct coordinates *)malloc(sizeof(struct coordinates)*Molecule.num_atoms);

   g = chkpt_rd_geom();
   Z = chkpt_rd_zvals();

   /*--- move it into the appropriate struct form ---*/
   for (i=0; i<Molecule.num_atoms; i++){
      Molecule.centers[i].x = g[i][0];
      Molecule.centers[i].y = g[i][1];
      Molecule.centers[i].z = g[i][2];
      Molecule.centers[i].Z_nuc = Z[i];
   }
   free_block(g);
   free(Z);

   return;
}


/*!---------------------------------------
  enuc computes nuclear repulsion energy
 ---------------------------------------*/

void compute_enuc()
{
  int i, j;
  double Z1Z2, r2, oor;
  double E = 0.0;

  if(Molecule.num_atoms > 1)
    for(i=1; i<Molecule.num_atoms; i++)
      for(j=0; j<i; j++){
	r2 = 0.0;
	r2 += (Molecule.centers[i].x-Molecule.centers[j].x)*
	     (Molecule.centers[i].x-Molecule.centers[j].x);
	r2 += (Molecule.centers[i].y-Molecule.centers[j].y)*
	     (Molecule.centers[i].y-Molecule.centers[j].y);
	r2 += (Molecule.centers[i].z-Molecule.centers[j].z)*
	     (Molecule.centers[i].z-Molecule.centers[j].z);
	oor = 1.0/sqrt(r2);
        Z1Z2 = Molecule.centers[i].Z_nuc*Molecule.centers[j].Z_nuc;
#ifdef __IBMCPP__
        if (isnan(oor) || isinf(oor)) {
#else
#ifdef HAVE_FUNC_ISINF
        if (std::isnan(oor) || std::isinf(oor)) {
#elif HAVE_FUNC_FINITE
        if (std::isnan(oor) || !std::finite(oor)) {
#endif
#endif
          if (fabs(Z1Z2) != 0.0)
            throw std::domain_error("compute_enuc -- charges too close to each other");
        }
        else
	  E += Z1Z2 * oor;
      }

  // include the contribution from the external electric field here
  // it doesn't fit anywhere else
  if (UserOptions.E_given) {
    // E field is given in the frame specified by EFIELD_FRAME
    // if necessary, rotate E to the canonical frame, in which all integrals are computed
    double* Ef = new double[3];  for(int i=0; i<3; ++i) Ef[i] = 0.0;
    switch(UserOptions.E_frame) {
      case reference:
      {
        double** rref = chkpt_rd_rref();
        for(int i=0; i<3; ++i)
          for(int j=0; j<3; ++j)
            Ef[i] += rref[i][j] * UserOptions.E[j];
        Chkpt::free(rref);
      }
      break;

      case canonical:
        for(int i=0; i<3; ++i) Ef[i] = UserOptions.E[i];
        break;

      default:
        throw std::runtime_error("This value for UserOptions.E_frame not supported. See documentation for keyword EFIELD_FRAME.");
    }
    double E_efield = 0.0;
    for(i=0; i<Molecule.num_atoms; i++) {
      E_efield -= Molecule.centers[i].Z_nuc * Molecule.centers[i].x * Ef[0];
      E_efield -= Molecule.centers[i].Z_nuc * Molecule.centers[i].y * Ef[1];
      E_efield -= Molecule.centers[i].Z_nuc * Molecule.centers[i].z * Ef[2];
    }
    delete[] Ef;

    fprintf(outfile,"\n    EFIELD(input file) = ( ");
    for(int i=0; i<3; ++i) fprintf(outfile," %12.9lf ",UserOptions.E[i]);
    fprintf(outfile," )\n");
    if (UserOptions.E_frame != canonical) {
      fprintf(outfile,"    EFIELD(%s frame) = ( ", (UserOptions.E_frame==canonical ? "canonical" : "reference"));
      for(int i=0; i<3; ++i) fprintf(outfile," %12.9lf ",Ef[i]);
      fprintf(outfile," )\n");
    }

    fprintf(outfile,"    Electric field contribution added to the nuclear repulsion energy:\n");
    fprintf(outfile,"      dE = %20.15lf\n", E_efield);

    E += E_efield;
  }

  Molecule.Enuc = E;

  return;
}

};};
