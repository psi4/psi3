/** FREQ_GRAD_CART computes frequencies from gradients and cartesian disps */

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
  #include <stdio.h>
  #include <libchkpt/chkpt.h>
  #include <stdlib.h>
  #include <string.h>
  #include <ctype.h>
  #include <libciomr/libciomr.h>
  #include <libqt/qt.h>
  #include <libipv1/ip_lib.h>
  #include <physconst.h>
  #include <libpsio/psio.h>
  #include <psifiles.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

extern double **compute_B(internals &simples, salc_set &salcs);
extern double *compute_q(internals &simples, salc_set &symm);
extern double **compute_G(double **B, int num_intcos, cartesians &carts);
extern int get_irrep_xyz(double **cartrep, int xyz);
void sort_evals_all(int nsalc_all, double *evals_all, int *evals_all_irrep);

/**** FREQ_GRAD_CART compute frequencies from gradients for cartesian
  displacements ****/

void freq_grad_cart(cartesians &carts) {
  int i,j,k,l,a,b, ii, cnt, dim, natom, xyz, cnt_eval = -1, *evals_all_irrep,op_disp;
  int xyzA, xyzB, atomA, atomB, row, col,h,nirreps,xyz_irr,I,cnt_all,match,sign;
  int nsalcs,ncoord,start_disp, start_salc, atom, atom2, op, loner, natom_unique;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert, k_convert;
  double *f, *f_q, *temp_arr, *temp_arr2, tval, **geom2D, **B_inv;
  double **evects, *evals, **FG, tmp, **force_constants_x, **tmp_mat_inv;
  double *micro_e, **micro_geom, **disp_grad, *grad, **grads_adapted, **tmp_mat;
  double *evals_all, **cartrep, **disp_grad_all;

  int *nsalc, *ndisp, ndisp_all, nsalc_all, **ict;
  double ***salc, ***disp;

  nirreps = syminfo.nirreps;
  natom = carts.get_natom();
  masses = carts.get_mass();
  u = mass_mat(masses);
  ndisp = init_int_array(nirreps);
  nsalc = init_int_array(nirreps);

  chkpt_init(PSIO_OPEN_OLD);
  cartrep = chkpt_rd_cartrep();
  ict = chkpt_rd_ict();
  natom_unique = chkpt_rd_num_unique_atom();
  chkpt_close();

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

  fprintf(outfile,"nsalc_all: %d ndisp_all%d\n", nsalc_all, ndisp_all);
  fprintf(outfile,"nsalc:"); for (h=0; h<nirreps; ++h) fprintf(outfile,"%d",nsalc[h]);
  fprintf(outfile,"ndisp:"); for (h=0; h<nirreps; ++h) fprintf(outfile,"%d",ndisp[h]);

  disp_grad = block_matrix(ndisp_all,3*natom);
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
    (char *) &(disp_grad[0][0]), ndisp_all*3*natom*sizeof(double));
  fprintf(outfile,"Gradients of displaced geometries\n");
  print_mat(disp_grad,ndisp_all,3*natom,outfile);

  /**** construct full gradient list using unique ones ****/
  disp_grad_all = block_matrix(2*nsalc_all, 3*natom);
  /* just copy symmetric gradients */
  for (j=0; j<ndisp[0]; ++j) {
    for (I=0; I<3*natom; ++I)
      disp_grad_all[j][I] = disp_grad[j][I];
  }

  /* add gradients to list for asymmetric displacements */
  cnt_all = cnt = ndisp[0] - 1;
  for (h=1; h<nirreps; ++h) {

    if (!ndisp[h]) continue;

    /* copy computed displacement in */
    for (j=0; j<ndisp[h]; ++j) {
      ++cnt;
      ++cnt_all;
      for (I=0; I<3*natom; ++I)
        disp_grad_all[cnt_all][I] = disp_grad[cnt][I];

    /* determine which operation takes minus displacement into plus displacement */
      for (op_disp=0; op_disp<syminfo.nirreps; ++op_disp) {
        if ( syminfo.ct[h][op_disp] == -1 )
          break;
        /*
        match = 1;
        for (atom=0; atom<natom_unique; ++atom) {
          atom2 = ict[op_disp][atom] - 1;
          for (xyz=0; xyz<3; ++xyz)
            if ( cartrep[op_disp][3*xyz+xyz] * disp_grad[cnt][3*atom2+xyz]
                                     != -1.0 * disp_grad[cnt][3*atom +xyz])
              match = 0;
        }
        if (match) break;
        */
      }
      fprintf(outfile,"Operation that takes plus displacement %d to minus is %s.\n",
          cnt, syminfo.op_lbls[op_disp]);

      /* what is the parity of the irrep for this special operation ? */
      sign = syminfo.ct[h][op_disp];
      fprintf(outfile,"Parity of irrep %s for this operation is %d.\n",
         syminfo.clean_irrep_lbls[h], sign);

      ++cnt_all;
      for (atom=0; atom<natom; ++atom) {
        for (xyz=0; xyz<3; ++xyz) {
          atom2 = ict[op_disp][atom] - 1;
          disp_grad_all[cnt_all][3*atom2+xyz] = /*sign **/ disp_grad[cnt][3*atom+xyz] * cartrep[op_disp][3*xyz+xyz];
        }
      }
    }
    ndisp[h] += ndisp[h];
  }
  free_block(disp_grad);

  /* fix number of displacements to match full list of gradients */
  ndisp_all=0;
  for (h=0; h<nirreps; ++h) {
    ndisp_all += ndisp[h];
  }
  fprintf(outfile,"total constructed ndisp_all: %d\n",ndisp_all);

  B = block_matrix(nsalc_all,3*natom);
  psio_read_entry(PSIF_OPTKING, "OPT: Adapted cartesians",
    (char *) &(B[0][0]), nsalc_all*3*natom*sizeof(double));
  fprintf(outfile,"B matrix\n");
  print_mat(B,nsalc_all,3*natom,outfile);

  evals_all = init_array(nsalc_all);
  evals_all_irrep = init_int_array(nsalc_all);

  /* not used 
  micro_e = init_array(ndisps);
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
    (char *) &(micro_e[0]), ndisps*sizeof(double));
  micro_geom = block_matrix(ndisps,3*natom);
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
    (char *) &(micro_geom[0][0]), ndisps*3*natom*sizeof(double));
  */

  /**** loop over irreps of salcs and vibrations ****/
  start_salc = 0;
  start_disp = 0;
  for (h=0; h<nirreps; ++h) {

    if (!nsalc[h]) continue;

    /* construct B_inv = (BuBt)-1 B u = G-1 B u */
    G = compute_G(&(B[start_salc]),nsalc[h],carts);
    G_inv = symm_matrix_invert(G,nsalc[h],1,0);
    free_block(G);

    /* mass-weight gradients; g_xm = (1/sqrt(m)) * g_x */
    /* works for H2 what is more general? */
    for (k=0; k<ndisp[h]; ++k)
      for (j=0;j<3*natom;++j)
        disp_grad_all[k+start_disp][j] /= sqrt(masses[j]);
    fprintf(outfile,"Displaced gradients/sqrt(masses[j])\n");
    print_mat(&(disp_grad_all[start_disp]),ndisp[h],3*natom,outfile);
    fflush(outfile);

    // compute forces in internal coordinates, f_q = G_inv B u f_x
    f_q = init_array(nsalc[h]);
    temp_arr = init_array(nsalc[h]);
    temp_arr2 = init_array(3*natom);
    grads_adapted = block_matrix(ndisp[h],3*natom);
  
    fprintf(outfile,"Gradients recomputed from internal coordinate gradients\n");
    for (k=0; k<ndisp[h]; ++k) {

      fprintf(outfile,"ndisp[h]: %d, start_disp: %d\n", ndisp[h], start_disp);
      mmult(u,0,&(disp_grad_all[k+start_disp]),1,&temp_arr2,1,3*natom,3*natom,1,0);
      mmult(&(B[start_salc]),0,&temp_arr2,1,&temp_arr,1,nsalc[h],3*natom,1,0);
      mmult(G_inv,0,&temp_arr,1,&f_q,1,nsalc[h],nsalc[h],1,0);
  
    for (j=0;j<3*natom;++j)
      grads_adapted[k][j] = f_q[j];

      /* test by transforming f_q back to cartesian forces and compare */
      mmult(B,1,&(grads_adapted[k]),1,&temp_arr2,1,3*natom,nsalc[h],1,0);
      print_mat2(&temp_arr2, 1, 3*natom, outfile);
    }

    free(f_q);
    free(temp_arr);
    free(temp_arr2);
    free_block(G_inv);

    /* fprintf(outfile,"Adapted gradients\n");
    print_mat(grads_adapted,ndisps,ncoord,outfile);
    fflush(outfile); */

    fprintf(outfile," ** Using %d-point formula.\n",optinfo.points_freq);
    force_constants = block_matrix(nsalc[h],nsalc[h]);

    /** construct force constant matrix **/
    for (i=0; i<nsalc[h]; ++i) {
      for (j=0; j<nsalc[h]; ++j) {
        force_constants[i][j] = 
          (grads_adapted[2*i+1][j] - grads_adapted[2*i][j])
            / (2.0 * optinfo.disp_size);
      }
    }
    fprintf(outfile,"Force Constants\n");
    print_mat(force_constants, nsalc[h], nsalc[h], outfile);
    fflush(outfile);

    start_salc += nsalc[h];
    start_disp += ndisp[h];

    dim = nsalc[h];

    /* masses = carts.get_mass();
    for (i=0;i<3*natom;++i)
      for (j=0;j<3*natom;++j) {
        force_constants[i][j] *= 1.0 / sqrt(masses[i] * masses[j]);
      }
    fprintf(outfile,"Mass Weighted Force Constants 3N by 3N\n");
    print_mat(force_constants, dim, dim, outfile);
    fflush(outfile); */

    /** Find eigenvalues of force constant matrix **/
    evals  = init_array(dim);
    evects = block_matrix(dim, dim);
    dgeev_optking(dim, force_constants, evals, evects);
    free_block(force_constants);
    sort(evals, evects, dim);
    free_block(evects);

    for (i=0; i<dim; ++i) {
      ++cnt_eval;
      evals_all[cnt_eval] = evals[i];
      evals_all_irrep[cnt_eval] = h;
    } 
    free(evals);
  }
  free(masses);
  free_block(u);
  free_block(disp_grad_all);

  sort_evals_all(nsalc_all,evals_all, evals_all_irrep);

  /* convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  /* v = 1/(2 pi c) sqrt( eval ) */
  fprintf(outfile, "\n\t        Harmonic Frequencies  \n");
  fprintf(outfile,   "\t              (cm-1)          \n");
  fprintf(outfile,   "\t------------------------------\n");
  k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for(i=nsalc_all-1; i>-1; --i) {
    if(evals_all[i] < 0.0)
      fprintf(outfile, "\t %3d %5s %10.3fi\n", i, syminfo.irrep_lbls[evals_all_irrep[i]],
          cm_convert * sqrt(-k_convert * evals_all[i]));
    else
      fprintf(outfile, "\t %3d %5s %10.3f\n", i, syminfo.irrep_lbls[evals_all_irrep[i]],
          cm_convert * sqrt(k_convert * evals_all[i]));
    }
  fprintf(outfile,   "\t----------------------------\n");
  fflush(outfile);
  free(evals_all);
  free(evals_all_irrep);

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();
  return;
}

void sort_evals_all(int nsalc_all, double *evals_all, int *evals_all_irrep) {
  int i,j,tmp_irrep,min_index;
  double min,tmp_eval;

  for (j=0; j<nsalc_all; ++j) {

    min = 1.0E6;
    for (i=j; i<nsalc_all; ++i) {
      if (evals_all[i] < min) {
        min = evals_all[i];
        min_index = i;
      }
    }
    tmp_eval = evals_all[j];
    evals_all[j]=evals_all[min_index];
    evals_all[min_index] = tmp_eval;

    tmp_irrep = evals_all_irrep[j];
    evals_all_irrep[j] = evals_all_irrep[min_index];
    evals_all_irrep[min_index] = tmp_irrep;
  }
  return;
}
