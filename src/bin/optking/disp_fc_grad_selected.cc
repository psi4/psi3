/*! \file disp_fc_grad_selected.cc
    \ingroup OPTKING
    \brief makes +/- displacements suitable for determining selected force constants from gradients
*/

#include <cmath>
#include <cstdio>
#include <libchkpt/chkpt.h>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

namespace psi { namespace optking {

extern double *compute_q(internals &simples, salc_set &symm);
extern int new_geom(cartesians &carts, internals &simples, salc_set &symm, 
    double *dq, int print_to_geom_file, int restart_geom_file,
    char *disp_label, int disp_num, int last_disp, double *return_geom);
extern int make_disp_irrep(cartesians &carts, internals &simples, salc_set &symm);

// only the symmetric salcs are passed in
// make displacements for selected coordinates

int disp_fc_grad_selected(cartesians &carts, internals &simples, salc_set &symm) 
{
  int i, j, errcod, ndisps, *irrep_per_disp, cnt, success;
  double *fgeom, energy, **micro_geoms, **displacements;
  //int simple_id, sub_index, sub_index2, intco_type;
  char disp_label[MAX_LINELENGTH];

  int ncoord = 0;        // # of coordinates to displace
  int nsymm = symm.get_num(); // # of symmetric coordinates
  int dim_carts = 3*carts.get_natom();
  int *coord2salc; // list of absolute salc indices for coordinates to displace

  errcod = ip_count("SELECTED_FC", &ncoord, 0);
  if (errcod != IPE_OK) throw("could not read selected_fc\n");
  coord2salc = new int [ncoord];
  for (i=0; i<ncoord; ++i) {
    errcod = ip_data("SELECTED_FC", "%d", &j, 1, i);
    coord2salc[i] = j-1;
    if (errcod != IPE_OK) throw("could not read selected_fc\n");
  }

  // assume all coordinates are totally symmetric and using 3-point formula
  ndisps = 2*ncoord;
  fprintf(outfile,"\nDoing %d displacements for %d coordinates.\n", ndisps, ncoord);
  fprintf(outfile,"Doing displacements for coordinates: ");
  for (i=0; i<ncoord; ++i) fprintf(outfile," %d", coord2salc[i]+1);

  // save reference geometry and energy
  fgeom = carts.get_coord();
  energy = carts.get_energy();
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry", (char *) fgeom, dim_carts*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy", (char *) &(energy), sizeof(double));
  close_PSIF();

/* automate interfragment coordinates later?
  // find number of single uncombined interfragment-type coodinates
  int ncoord=0;
  for (i=0; i<nsymm; ++i) {
    simple_id = symm.get_simple(i,0);
    simples.locate_id(simple_id,&intco_type,&sub_index,&sub_index2);
    if ( (symm.get_length(i) != 1) || (intco_type != FRAG_TYPE) ) continue;
    ++ncoord;
  }
  cnt = 0;
  for (i=0; i<ncoord; ++i) {
    simple_id = symm.get_simple(i,0);
    simples.locate_id(simple_id,&intco_type,&sub_index,&sub_index2);
    if ( (symm.get_length(i) != 1) || (intco_type != FRAG_TYPE) ) continue;
    coord2salc[cnt] = i;
  }
*/

  displacements = block_matrix(ndisps, nsymm);
  for (i=0; i<ncoord; ++i) {
    displacements[2*i  ][coord2salc[i]] = -1.0 * optinfo.disp_size;
    displacements[2*i+1][coord2salc[i]] = 1.0 * optinfo.disp_size;
  }
  if (optinfo.print_hessian) {
    fprintf(outfile,"\nDisplacement Matrix\n");
    print_mat5(displacements, ndisps, nsymm, outfile);
  }

  /*** generate and store Micro_iteration cartesian geometries ***/
  micro_geoms = block_matrix(ndisps, dim_carts);
  for (i=0;i<ndisps;++i)  {
    sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
    success = new_geom(carts,simples,symm,displacements[i],0,
        0, disp_label, i, 0, micro_geoms[i]);
    if (!success) {
      fprintf(outfile,"Unable to generate displaced geometry.\n");
      exit_io();
      exit(PSI_RETURN_FAILURE);
    }
  }
  free_block(displacements);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geoms[0][0]), ndisps*dim_carts*sizeof(double));

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisps), sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced coords to Salc Index",
      (char *) coord2salc, ncoord * sizeof(int));

  close_PSIF();
  free_block(micro_geoms);
  free(coord2salc);

  // write zeroes for initial energy and gradients of displacements
  double *disp_e, *disp_grad;

  disp_e = new double[ndisps];
  disp_grad = new double [ndisps*dim_carts*sizeof(double)];

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(disp_grad[0]), ndisps*dim_carts*sizeof(double));

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), ndisps*sizeof(double));

  // assume all displacements are totally symmetric
  irrep_per_disp = init_int_array(ndisps);
  psio_write_entry(PSIF_OPTKING, "OPT: Irrep per disp",
    (char *) &(irrep_per_disp[0]), ndisps*sizeof(int));
  free(irrep_per_disp);

  close_PSIF();

  // Reset microiteration counter
  optinfo.micro_iteration = 0;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Micro_iteration",
      (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();

  delete [] disp_e;
  delete [] disp_grad;
  return(ndisps);
}

}} /* namespace psi::optking */

