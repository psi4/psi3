/*! \file
    \ingroup OPTKING
    \brief grad_energy_cart(): compute gradients from displacements
   in cartesian coordinates - displacements are done by disp_grad_energy_cart.cc
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libqt/qt.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>

namespace psi { namespace optking {

void grad_energy_cart(cartesians & carts) {
  int i, j, k, l, dim, natom, cnt_eval = -1, *evals_all_irrep;
  int h, nirreps, *nsalc, *ndisp, ndisp_all, nsalc_all, *start_irr;
  double **B, **force_constants, energy_ref, *energies, cm_convert, k_convert;
  double *f, tval, **evects, *evals, tmp, disp_size;
  double *disp_E, *evals_all, **cartrep, ***salc, ***disp, **normal;
  int print;
  char *line1;
  print = optinfo.print_cartesians;

  nirreps = syminfo.nirreps;
  natom = optinfo.natom;
  disp_size = optinfo.disp_size;

  chkpt_init(PSIO_OPEN_OLD);
  cartrep = chkpt_rd_cartrep();
  chkpt_close();

  ndisp = init_int_array(nirreps);
  nsalc = init_int_array(nirreps);

  /* Read in data from PSIF_OPTKING */
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisp_all), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Num. of coord.",
      (char *) &(nsalc_all), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Disp. per irrep",
      (char *) &(ndisp[0]), nirreps*sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Coord. per irrep",
      (char *) &(nsalc[0]), nirreps*sizeof(int));

  fprintf(outfile,"total nsalc: %d, total ndisp: %d\n", nsalc_all, ndisp_all);
  fprintf(outfile,"nsalc per irrep: "); for (h=0; h<nirreps; ++h) fprintf(outfile,"%d ",nsalc[h]);
  fprintf(outfile,"\n");
  fprintf(outfile,"ndisp per irrep: "); for (h=0; h<nirreps; ++h) fprintf(outfile,"%d ",ndisp[h]);
  fprintf(outfile,"\n\n");

  disp_E = init_array(ndisp_all);

  if (optinfo.energy_dat) { /* read energy.dat text file */
    fprintf(outfile,"Reading displaced energies from energy.dat.\n");
    fp_energy_dat = fopen("energy.dat", "r");
    if (fp_energy_dat == NULL) {
      fprintf(outfile,"energy.dat file not found.\n");
      exit(PSI_RETURN_FAILURE);
    }
    rewind (fp_energy_dat);
    line1 = new char[MAX_LINELENGTH+1];
    // ACS (11/06) Allow external program to be used to compute energies
    if(optinfo.external_energies){
      /* Read the first energy and dump it as the reference energy */
      double temp;
      fscanf(fp_energy_dat,"%lf",&temp);
      psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",(char *) &(temp), sizeof(double));
    }
    for (i=0; i<ndisp_all; i++) {
      fscanf(fp_energy_dat,"%lf",&(disp_E[i]));
    }

    fclose(fp_energy_dat);
    delete [] line1;
  }
  else {
    fprintf(outfile,"Reading displaced energies from PSIF_OPTKING.\n");
    psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
        (char *) &(disp_E[0]), ndisp_all*sizeof(double));
  }

  fprintf(outfile,"Energies of displaced geometries. Check for precision!\n");
  int cnt = -1;
  for (i=0;i<nsalc[0];++i) {
    fprintf(outfile,"Coordinate %d: ",i);
    for (j=0;j<optinfo.points-1;++j)
      fprintf(outfile,"%15.10lf ", disp_E[++cnt]);
    fprintf(outfile,"\n");
  }
  fflush(outfile);

  // Calculate forces in mass-weighted, symmetry-adapted cartesians
  double *f_q = new double[nsalc[0]];
  if (optinfo.points == 3) {
    for (i=0; i<nsalc[0]; ++i) {
      f_q[i] = (disp_E[2*i+1]-disp_E[2*i]) / (2.0 * optinfo.disp_size);
      f_q[i] = -1.0 * f_q[i] * _hartree2J * 1.0E18 ;
    }
  }
  else if (optinfo.points == 5) {
    for (i=0; i<nsalc[0]; ++i) {
      f_q[i] = ( disp_E[4*i]-8.0*disp_E[4*i+1]+8.0*disp_E[4*i+2]-disp_E[4*i+3])
                  / (12.0 * optinfo.disp_size);
      f_q[i] = -1.0 * f_q[i] * _hartree2J * 1.0E18 ;
    }
  }

  fprintf(outfile,"\nForces in symmetric cartesian coordinates:\n");
  for (i=0; i<nsalc[0]; ++i)
    fprintf(outfile,"%13.10lf\n",f_q[i]);

  B = init_matrix(nsalc_all,3*natom);

  double *geom = new double[3*natom];

  psio_read_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy_ref), sizeof(double));
  psio_read_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(geom[0]), 3*natom*sizeof(double));
  psio_read_entry(PSIF_OPTKING, "OPT: Adapted cartesians",
    (char *) &(B[0][0]), nsalc_all*3*natom*sizeof(double));
  close_PSIF();

  if (print) {
    fprintf(outfile,"Reference energy: %15.10lf\n",energy_ref);
    fprintf(outfile,"B matrix (adapted cartesians)\n");
    print_mat(B,nsalc[0],3*natom,outfile);
  }


// Transform forces back into ordinary cartesians.
  //simples.compute(geom);
  //simples.compute_s(geom);
  //B = compute_B(simples, symm);

  // B matrix (read above) is nsalc by 3*natom
  f = new double[3*natom];
  opt_mmult(B,1,&f_q,1,&f,1,3*natom, nsalc[0],1,0);
  free_matrix(B);

  // change forces to gradient for writing a file11 entry
  for(i=0;i<3*natom;++i)
    f[i] = -1.0 * f[i] / _hartree2J / 1.0E18 * _bohr2angstroms;

  // write out file11.dat
  char *disp_label = new char[MAX_LINELENGTH];

  opt_ffile(&fp_11, "file11.dat", 1);
  char *wfn,*dertype;
  int errcod = ip_string("WFN", &wfn, 0);
  if (errcod != IPE_OK)
    punt("Keyword WFN not found in input file");
  errcod = ip_string("DERTYPE", &dertype, 0);
  if (errcod != IPE_OK)
    dertype = strdup("NONE");
  chkpt_init(PSIO_OPEN_OLD);
  char* label = chkpt_rd_label();
  chkpt_close();
  sprintf(disp_label,"%-59.59s %-10.10s%-8.8s",label,wfn,dertype);
  free(label); free(wfn); free(dertype);
  carts.set_energy(energy_ref);
  carts.set_coord(geom);
  carts.set_grad(f);
  carts.print(11,fp_11,0,disp_label, 0);
  fclose(fp_11);

  // write out geometry, gradient and energy to chkpt file
  cnt = -1;
  double **geom2D = init_matrix(carts.get_natom(),3);
  for (i=0; i<carts.get_natom(); ++i)
    for (j=0; j<3; ++j)
      geom2D[i][j] = geom[++cnt];

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom2D);
  chkpt_wt_grad(f);
  chkpt_wt_etot(energy_ref);
  chkpt_close();
  free_matrix(geom2D);
  free_array(f);
  free_array(geom);

  return;
}

}}

