/*! \file
    \ingroup OPTKING
    \brief displaces along all coordinates + and - by disp_size
    need to save geometry to last step for gradients by energy
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
extern int new_geom(cartesians &carts, internals &simples, salc_set &all_salcs, 
    double *dq, int print_to_geom_file, int restart_geom_file,
    char *disp_label, int disp_num, int last_disp, double *return_geom);

/* MAKE_DISP_IRREP - make displacements for modes labelled with 
* symmetry IRREP (+ and - if symmetric; otherwise, just +) */

int make_disp_irrep(cartesians &carts, internals &simples, salc_set &all_salcs) 
{
  int i,j,a,b, cnt,dim, dim_carts, ndisps, nsalcs, *irrep_salcs, irrep;
  int *irrep_per_disp, success;
  double *fgeom, energy, **micro_geoms, **displacements;
  double *f, *q, tval;
  char *disp_label, **disp_irrep_lbls, *salc_lbl;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_natom();
  irrep_salcs = new int[all_salcs.get_num()];
  irrep = optinfo.irrep;

  /* count number of IRREP salcs */
  nsalcs = 0;
  cnt = 0;
  for (i=0; i<all_salcs.get_num(); ++i) {
    salc_lbl = all_salcs.get_label(i);
    if ( strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0) {
      ++nsalcs;
      irrep_salcs[cnt++] = i;
    }
  }
  fprintf(outfile,"Found %d internals of irrep %d \n", nsalcs, irrep+1);
//  for (i=0; i<nsalcs; ++i)
//    fprintf(outfile,"internal id %d\n",irrep_salcs[i]);

  if (nsalcs == 0) { 
    fprintf(outfile,"Produced 0 displacements.\n");
    return(0);
  }

  // write reference geometry (the one being replaced to PSIF)
  fgeom = carts.get_coord();
  energy = carts.get_energy();
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(fgeom[0]), dim_carts*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  /*** make list of internal displacements for micro_iterations ***/
  if (irrep == 0) {
    if (optinfo.points == 3) {
      ndisps = 2*nsalcs;
      displacements = block_matrix(ndisps,all_salcs.get_num());
      cnt = 0;
      for (i=0; i<nsalcs; ++i) {
        displacements[cnt++][irrep_salcs[i]] = -1.0 * optinfo.disp_size;
        displacements[cnt++][irrep_salcs[i]] =  1.0 * optinfo.disp_size;
      }
    }
    else if (optinfo.points == 5) {
      ndisps = 4*nsalcs;
      displacements = block_matrix(ndisps,all_salcs.get_num());
      cnt = 0;
      for (i=0; i<nsalcs; ++i) {
        displacements[cnt++][irrep_salcs[i]] = -2.0 * optinfo.disp_size;
        displacements[cnt++][irrep_salcs[i]] = -1.0 * optinfo.disp_size;
        displacements[cnt++][irrep_salcs[i]] =  1.0 * optinfo.disp_size;
        displacements[cnt++][irrep_salcs[i]] =  2.0 * optinfo.disp_size;
      }
    }
  }
  else { // non-symmetric irrep
    if (optinfo.points == 3) {
      ndisps = nsalcs;
      displacements = block_matrix(ndisps,all_salcs.get_num());
      cnt = 0;
      for (i=0; i<nsalcs; ++i)
        displacements[cnt++][irrep_salcs[i]] = -1.0 * optinfo.disp_size;
    }
    else if (optinfo.points == 5) {
      ndisps = 2*nsalcs;
      displacements = block_matrix(ndisps,all_salcs.get_num());
      cnt = 0;
      for (i=0; i<nsalcs; ++i) {
        displacements[cnt++][irrep_salcs[i]] = -2.0 * optinfo.disp_size;
        displacements[cnt++][irrep_salcs[i]] = -1.0 * optinfo.disp_size;
      }
    }
  }

  fprintf(outfile,"Generating a total of %d displacements ", ndisps);
  fprintf(outfile,"using %d-point formula for modes of irrep %d.\n",
      optinfo.points, irrep+1);

  fprintf(outfile,"\nDisplacement Matrix\n");
  print_mat5(displacements, ndisps, all_salcs.get_num(), outfile);

  /*** generate and store Micro_iteration cartesian geometries ***/
  micro_geoms = block_matrix(ndisps, dim_carts);
  if (optinfo.freeze_intrafragment) { // compute new geometry analytically
    int simple, intco_type, sub_index, sub_index2, J, simple_b, nf, disp, xyz, frag;
    int A_natom, B_natom, *A_atom, *B_atom, *frag_atom, *frag_allatom;
    double **geom_A, **geom_B, **geom_2D;
    double R_AB, theta_A, theta_B, tau, phi_A, phi_B;

    // make sure that all coordinates are simple and inter-fragment
    for (i=0;i<all_salcs.get_num();++i) {
      if (all_salcs.get_length(i) != 1) {
        printf("Only simple coordinates can be used with freeze_intrafragment\n");
        exit(PSI_RETURN_FAILURE);
      }
      simple = all_salcs.get_simple(i,0);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);
      if (intco_type != FRAG_TYPE) {
        printf("Only inter-fragment coordinates can be used with freeze_intrafragment\n");
        exit(PSI_RETURN_FAILURE);
      }
    }

    nf = optinfo.nfragment; // move to getoptinfo at some point
    frag_atom = (int *) malloc(nf*sizeof(int));
    frag_allatom = (int *) malloc(nf*sizeof(int));
    frag_atom[0] = frag_allatom[0] = 0;
    for (frag=1; frag<nf; ++frag) {
      frag_atom[frag] = frag_atom[frag-1] + optinfo.natom_per_fragment[frag-1];
      frag_allatom[frag] = frag_allatom[frag-1] + optinfo.natom_per_fragment[frag-1];
    }
  
    for (disp=0;disp<ndisps;++disp)  { // loop over displacements
      geom_2D = carts.get_coord_2d(); // in bohr

      //load first fragment
      a = optinfo.natom_per_fragment[0];
      for (i=0; i<a; ++i)
        for (xyz=0; xyz<3; ++xyz)
            micro_geoms[disp][3*i+xyz] = geom_2D[i][xyz];

      cnt = all_salcs.get_num();
      for (frag=1; frag<optinfo.nfragment; ++frag) { // loop over remaining fragments
        a = optinfo.natom_per_fragment[frag-1];
        b = optinfo.natom_per_fragment[frag];
        geom_A = block_matrix(a,3);
        geom_B = block_matrix(b,3);
        for (xyz=0;xyz<3;++xyz) {
          for (i=0;i<a;++i)
            geom_A[i][xyz] = geom_2D[frag_atom[frag-1]+i][xyz];
          for (i=0;i<b;++i)
            geom_B[i][xyz] = geom_2D[frag_atom[frag]+i][xyz];
        }

        q = compute_q(simples,all_salcs);

        R_AB = theta_A = theta_B = tau = phi_A = phi_B = 0.0;
        if (cnt > 0) R_AB    = q[6*(frag-1)+0] + displacements[disp][6*(frag-1)+0];
        --cnt;
        if (cnt > 0) theta_A = q[6*(frag-1)+1] + displacements[disp][6*(frag-1)+1];
        --cnt;
        if (cnt > 0) theta_B = q[6*(frag-1)+2] + displacements[disp][6*(frag-1)+2];
        --cnt;
        if (cnt > 0) tau     = q[6*(frag-1)+3] + displacements[disp][6*(frag-1)+3];
        --cnt;
        if (cnt > 0) phi_A   = q[6*(frag-1)+4] + displacements[disp][6*(frag-1)+4];
        --cnt;
        if (cnt > 0) phi_B   = q[6*(frag-1)+5] + displacements[disp][6*(frag-1)+5];
        --cnt;

        R_AB /= _bohr2angstroms;
        theta_A *= 180.0/_pi;
        theta_B *= 180.0/_pi;
        tau *= 180.0/_pi; 
        phi_A *= 180.0/_pi; 
        phi_B *= 180.0/_pi; 

        orient_fragment(a, b, optinfo.nref_per_fragment[frag-1], optinfo.nref_per_fragment[frag],
          geom_A, geom_B, optinfo.fragment_coeff[frag-1], optinfo.fragment_coeff[frag],
          R_AB, theta_A, theta_B, tau, phi_A, phi_B, outfile);

        for (i=0; i<b; ++i)
          for (xyz=0; xyz<3; ++xyz)
            micro_geoms[disp][3*(frag_atom[frag]+i)+xyz] = geom_B[i][xyz];

         free_block(geom_A);
         free_block(geom_B);
      }
      free_block(geom_2D);
    }
    free(frag_atom);
    free(frag_allatom);
  }
  else {  /* do iterative back-transformation */
    for (i=0;i<ndisps;++i)  {
      sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
      success = new_geom(carts,simples,all_salcs,displacements[i],0,
          0, disp_label, i, 0, micro_geoms[i]);
      if (!success) {
        fprintf(outfile,"Unable to generate displaced geometry.\n");
        fclose(outfile);
        exit(PSI_RETURN_FAILURE);
      }
    }
  }
  free_block(displacements);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geoms[0][0]), ndisps*dim_carts*sizeof(double));

  if(optinfo.external_energies){
    FILE *fp_dispcart = fopen("dispcart", "w");
    // ACS (11/06) Print the displaced and reference geometry for external programs
    // Use 1 based counting, and remember that the first "displacement" is the reference itself
    fprintf(fp_dispcart,"%d\n",1);
    for(int n = 0; n < carts.get_natom(); n++){
      fprintf(fp_dispcart,"%16.10lf %16.10lf %16.10lf\n",fgeom[3*n],fgeom[3*n+1],fgeom[3*n+2]);
    }

    for(int disp = 0; disp< ndisps; disp++){
      fprintf(fp_dispcart,"%d\n",disp+2);
      for(int n = 0; n < carts.get_natom(); n++){
        fprintf(fp_dispcart,"%16.10lf %16.10lf %16.10lf\n",micro_geoms[disp][3*n],micro_geoms[disp][3*n+1],micro_geoms[disp][3*n+2]);
      }
    }
    fclose(fp_dispcart);
  }

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Num. of disp.",
      (char *) &(ndisps), sizeof(int));

  free_block(micro_geoms);

  // make room for storage of energy and gradients of displacements
  double *disp_e, *disp_grad;

  disp_e = new double[ndisps];
  disp_grad = new double [ndisps*3*carts.get_natom()*sizeof(double)];
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(disp_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), ndisps*sizeof(double));

  irrep_per_disp = init_int_array(ndisps);
  for (i=0; i<ndisps; ++i) irrep_per_disp[i] = irrep;
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

  fprintf(outfile,"Produced a total of %d displacements.\n",ndisps);

  delete [] disp_label;
  delete [] irrep_salcs;
  delete [] disp_e;
  delete [] disp_grad;
  return(ndisps);
}


/* MAKE_DISP_NOSYMM generate displacements - do positive and
* negative displacements along all coordinates ignorning symmetry */

int make_disp_nosymm(cartesians &carts, internals &simples, salc_set &all_salcs) 
{
  int i,j,a,b, dim, dim_carts, ndisps, nsalcs, *irrep_per_disp, cnt, success;
  double *fgeom, energy, **micro_geoms, **displacements;
  double *f, *q, tval;
  char *disp_label, **disp_irrep_lbls, *lbl;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_natom();

  nsalcs = all_salcs.get_num();

  if (nsalcs == 0)
    punt("There are no appropriate SALCs present to displace along.\n");

  // write reference geometry (the one being replaced to PSIF)
  fgeom = carts.get_coord();
  energy = carts.get_energy();
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(fgeom[0]), dim_carts*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  /*** make list of internal displacements for micro_iterations ***/
  if (optinfo.points == 3) {
    ndisps = 2*nsalcs;
    displacements = block_matrix(ndisps, nsalcs);
    for (i=0;i<all_salcs.get_num();++i) {
      displacements[2*i][i] = -1.0 * optinfo.disp_size;
      displacements[2*i+1][i] = 1.0 * optinfo.disp_size;
    }
  }
  else if (optinfo.points == 5) {
    ndisps = 4*nsalcs;
    displacements = block_matrix(ndisps, nsalcs);
    for (i=0;i<all_salcs.get_num();++i) {
      displacements[4*i+0][i] = -2.0 * optinfo.disp_size;
      displacements[4*i+1][i] = -1.0 * optinfo.disp_size;
      displacements[4*i+2][i] =  1.0 * optinfo.disp_size;
      displacements[4*i+3][i] =  2.0 * optinfo.disp_size;
    }
  }

  fprintf(outfile,"Generating a total of %d displacements ", ndisps);
  fprintf(outfile,"using %d-point formula for all modes.\n", optinfo.points);

  fprintf(outfile,"\nDisplacement Matrix\n");
  print_mat5(displacements, ndisps, nsalcs, outfile);

  /*** generate and store Micro_iteration cartesian geometries ***/
  micro_geoms = block_matrix(ndisps, dim_carts);
  for (i=0;i<ndisps;++i)  {
    sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
    success = new_geom(carts,simples,all_salcs,displacements[i],0,
        0, disp_label, i, 0, micro_geoms[i]);
    if (!success) {
      fprintf(outfile,"Unable to generate displaced geometry.\n");
      fclose(outfile);
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

  close_PSIF();
  free_block(micro_geoms);

  // write zeroes for initial energy and gradients of displacements
  double *disp_e, *disp_grad;

  disp_e = new double[ndisps];
  disp_grad = new double [ndisps*3*carts.get_natom()*sizeof(double)];

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(disp_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), ndisps*sizeof(double));

  irrep_per_disp = init_int_array(ndisps);

  cnt = -1;
  for (i=0; i<ndisps; ++i) {
    lbl = all_salcs.get_label(i); /* returns pointer */
    for (j=0; j<syminfo.nirreps; ++j) {
      if ( strcmp( lbl, syminfo.irrep_lbls[j]) == 0) {
        irrep_per_disp[++cnt] = j;
        irrep_per_disp[++cnt] = j;
      }
    }
    // lbl pointer is freed when all_salcs gets deleted
  } 

  fprintf(outfile,"Irrep per displacement:\n");
  for (i=0;i<ndisps;++i)
    fprintf(outfile,"%3d",irrep_per_disp[i]);
  fprintf(outfile,"\n");

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


/*
void disp_docc(char **salc_lbl, int disp_nirrep, int *disp_clsdpi,
    int *disp_openpi, int *disp_frdocc, int *disp_fruocc) {
  int irrep, disp_irrep;
  char *ptgrp;

  ptgrp = syminfo.symmetry;

  for (irrep=0; irrep<syminfo.nirreps; ++irrep) {
    if (strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0) break;
  }
  fprintf(outfile,"Irrep of displacement %s or %d ", salc_lbl, irrep);

  static int nirrep_C1[1]  = {1};
  static int nirrep_CS[2]  = {2, 1};
  static int nirrep_C2V[4] = {4, 2, 2, 2};
  static int nirrep_C2H[4] = {4, 2, 2, 2};
  static int nirrep_D2H[8] = {8, 4, 4, 4, 4, 4, 4, 4};

  static int corr_table_C2V[4][2][4] = {
    { { 0, 0, 0, 0}, { 0, 0, 0, 0} }, // don't use for A1
    { { 1, 1, 0, 0}, { 0, 0, 1, 1} },
    { { 1, 0, 1, 0}, { 0, 1, 0, 1} },
    { { 1, 0, 0, 1}, { 0, 1, 1, 0} }
  }
  static int corr_table_CS[2][2][2] = {
    { { 0, 0} }, // don't use for A1
    { { 1, 1} }
  }
  static int corr_table_C2H[4][2][4] = {
    { { 0, 0, 0, 0}, { 0, 0, 0, 0} }, // don't use for A1
    { { 1, 1, 0, 0}, { 0, 0, 1, 1} },
    { { 1, 0, 1, 0}, { 0, 1, 0, 1} },
    { { 1, 0, 0, 1}, { 0, 1, 1, 0} }
  }
  static int corr_table_D2H[4][2][4] = {
    { { 0, 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0, 0} },
    { { 1, 1, 0, 0, 0, 0, 0, 0}, { 0, 0, 1, 1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0, 0} },
    { { 1, 0, 1, 0, 0, 0, 0, 0}, { 0, 1, 0, 1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0, 0} },
    { { 1, 0, 0, 1, 0, 0, 0, 0}, { 0, 1, 1, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0, 0} }
  }

  if (strcmp(ptgrp,"C1 ") == 0) {
    disp_nirrep = 1;
  }
  else if ((strcmp(ptgrp,"CS ") == 0) || (strcmp(ptgrp,"CI ") == 0) ||
      (strcmp(ptgrp,"C2 ") == 0)) {
    disp_nirrep = disp_nirrep_CS[irrep];
  }
  else if ((strcmp(ptgrp,"D2 ") == 0) || (strcmp(ptgrp,"C2V") == 0)) {
    disp_nirrep = disp_nirrep_C2V[irrep];
  }
  else if (strcmp(ptgrp,"C2H") == 0) {
    disp_nirrep = disp_nirrep_C2H[irrep];
  }
  else if (strcmp(ptgrp,"D2H") == 0) {
    disp_nirrep = disp_nirrep_D2H[irrep];
  }

  disp_clsdpi = new int[disp_nirrep];
  disp_openpi = new int[disp_nirrep];
  disp_frdocc = new int[disp_nirrep];
  disp_fruocc = new int[disp_nirrep];

  fprintf(outfile,"disp_nirrep: %d\n", disp_nirrep);
  return;
}
*/
}} /* namespace psi::optking */

