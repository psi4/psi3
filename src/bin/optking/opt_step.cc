/*** OPT_STEP.CC takes geometry steps using gradients -- optking's default operation ***/ 

#include <cmath>
extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <libciomr/libciomr.h>
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

extern double *compute_q(internals &simples, salc_set &symm);
extern double **compute_B(internals &simples, salc_set &symm);
extern double **compute_G(double **B, int num_intcos, cartesians &carts);
extern void new_geom(cartesians &carts, internals &simples, salc_set &symm, double *dq,
    int print_to_geom_file, int restart_geom_file, 
    char *disp_label, int disp_num, int last_disp, double *return_geom);
extern void empirical_H(internals &simples, salc_set &symm, cartesians &carts);
extern double **compute_H(salc_set &symm,double *q, double *f_q, double **P);
void fconst_init(cartesians &carts, internals &simples, salc_set &symm);
extern void compute_zmat(cartesians &carts, int *unique_zvars);
extern void print_zmat(FILE *outfile, int *unique_zvars);

int opt_step(cartesians &carts, internals &simples, salc_set &symm) {

  int i,j,a,b, dim, nallatom, dim_carts;
  double **B, **G, **G2, **G_inv, **H_inv, **temp_mat, **u, **P;
  double *temp_arr2, *temp_arr, *masses, *geom, *forces;
  double *f, *f_q, *dq, *q, tval, tval2, scale, temp, *djunk;
  char *disp_label, *wfn;

  // dim_carts = carts.get_natom()*3;
  nallatom = optinfo.nallatom;
  dim_carts = 3*nallatom;
  djunk = new double[dim_carts];
  disp_label = new char[MAX_LINELENGTH];

  if (symm.get_num() == 0)
    punt("No symmetric internal coordinates to optimize.\n");

  ip_string("WFN", &(wfn),0);
  fprintf(outfile,"\nCurrent %s energy before step %20.10lf\n",
      wfn, carts.get_energy());
  fprintf(outfile,"\nTaking geometry step number %d\n",optinfo.iteration+1);

  //*** Build transformation matrices
  dq = init_array(symm.get_num());
  q = compute_q(simples,symm);

  // build G = BuB^t
  B = compute_B(simples,symm);
  G = compute_G(B,symm.get_num(),carts);

  // compute G_inv
  fprintf(outfile,"\nBuB^t ");
  G_inv = symm_matrix_invert(G,symm.get_num(),1,optinfo.redundant);

  // setup the masses matrix, u
  masses = carts.get_fmass();
  u = mass_mat(masses);
  free(masses);

  // get forces array in cartesian coordinates, f, (in aJ/Ang)
  f = carts.get_fforces();

  //    fprintf(outfile,"cartesian forces in aJ/Ang");
  //    print_mat2(&f, 1, 3*carts.get_natom(), outfile);

  // compute forces in internal coordinates, f_q = G_inv B u f
  f_q = init_array(symm.get_num());
  temp_arr = init_array(symm.get_num());
  temp_arr2 = init_array(dim_carts);

  mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
  mmult(B,0,&temp_arr2,1,&temp_arr,1,symm.get_num(),dim_carts,1,0);
  mmult(G_inv,0,&temp_arr,1,&f_q,1,symm.get_num(),symm.get_num(),1,0);

  //fprintf(outfile,"internal forces in aJ/A");
  //print_mat2(&f_q, 1, symm.get_num(), outfile);
  // for (i=0;i<1;++i)
  //    f_q[i] = f_q[i] * sqrt(2) * _hartree2J * 1.0E18 / _bohr2angstroms;
  //    f_q[1] = f_q[1] * (_hartree2J * 1.0E18);

  free(f);
  free(temp_arr);
  free(temp_arr2);
  free_block(u);

  // Setup projection matrix P = G * G-
  // for inversion of Hessian with redundant internals
  P = block_matrix(symm.get_num(),symm.get_num());
  mmult(G,0,G_inv,0,P,0,symm.get_num(),symm.get_num(),symm.get_num(),0); 
  free_block(G);
  free_block(G_inv);

  // make sure some force constant are in PSIF
  fconst_init(carts, simples, symm);

  // Read in Hessian and update it if necessary from opt.aux
  H_inv = compute_H(symm,q,f_q,P);
  free_block(P);

  // Write Values and Forces of internals to opt.aux for later
  if (optinfo.bfgs) {
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Previous Internal Values",
        (char *) &(q[0]), symm.get_num()* sizeof(double));
    psio_write_entry(PSIF_OPTKING, "Previous Internal Forces",
        (char *) &(f_q[0]), symm.get_num()* sizeof(double));
    close_PSIF();
  }

  // Computing internal coordinate displacements H_inv f = dq, and new q
  mmult(H_inv,0,&f_q,1,&dq,1,symm.get_num(),symm.get_num(),1,0);
  free_block(H_inv);

  /* determine scale factor needed to keep step less than 10% of q if q big
     or less than 0.1 if q not so big, hopefully better solution coming soon */
  scale = 1;
  temp = 1;
  double cut = STEP_LIMIT / STEP_PERCENT;
  for (i=0;i<symm.get_num();++i) {
    if (fabs(dq[i]) > STEP_LIMIT) {
      if (fabs(dq[i]) > STEP_PERCENT*fabs(q[i]) && fabs(q[i]) > cut) { 
        temp = STEP_PERCENT*fabs(q[i])/fabs(dq[i]);
      }
      else if (fabs(q[i]) < fabs(dq[i]) || fabs(dq[i]) < cut) {
        temp = STEP_LIMIT / fabs(dq[i]);
      }
    }
    if (temp < scale){
      scale = temp;
    }
  }
  fprintf(outfile,"\nScaling displacements by %lf\n",scale); 
  for (i=0;i<symm.get_num();++i) {
    dq[i] = dq[i] * scale;   
    q[i] += dq[i];
  }

  // Print step summary to output.dat
  fprintf(outfile,
      "\nInternal Coordinate Update in Ang or Rad, aJ/Ang or aJ/Rad\n");
  fprintf(outfile,
      "         Value          Force          Displacement   New Value\n");
  for (i=0;i<symm.get_num();++i)
//    fprintf(outfile,"%2d%15.6lf%15.6lf%15.6lf%15.6lf\n",
    fprintf(outfile,"%2d%15.10lf%15.10lf%15.10lf%15.10lf\n",
        i+1,  q[i]-dq[i],    f_q[i],           dq[i],         q[i]);

  // Compute Max and RMS force, and see if geometry is optimized
  tval = 0.0;
  tval2 = fabs(f_q[0]);
  for (i=0;i<symm.get_num();++i) {
    tval += SQR(f_q[i]);
    if (fabs(f_q[i]) > tval2) tval2 = fabs(f_q[i]);
  }
  tval = tval/((double) symm.get_num());
  tval = sqrt(tval);
  fprintf(outfile,"   MAX force: %15.10lf   RMS force: %15.10lf\n",tval2,tval);
  if (tval2 < optinfo.conv) {
    fprintf(outfile,"\nMAX force is < %5.1e.  Optimization is complete.\n",
            optinfo.conv);
    fprintf(outfile,"Final %s energy is %15.10lf\n", wfn, carts.get_energy());
    fprintf(stderr,"\n  OPTKING:  optimization is complete\n");
    fprintf(outfile,"The Optimized geometry in a.u.\n");
    if (optinfo.zmat) {
      int *unique_zvars;
      unique_zvars = (int *) malloc(MAX_ZVARS*sizeof(int));
      compute_zmat(carts, unique_zvars);
      print_zmat(outfile, unique_zvars);
      free(unique_zvars);
      fprintf(outfile,"\n");
    }
    carts.print(12,outfile,0,disp_label,0);
    // fprintf(stderr,"\n  Returning code %d\n", PSI_RETURN_ENDLOOP);
    return(PSI_RETURN_ENDLOOP);
  }
  free(wfn);
  free(f_q);

  strcpy(disp_label,"New Cartesian Geometry in a.u.");
  new_geom(carts,simples,symm,dq,32,0,disp_label,0,0,djunk);
  free(q); free(dq);
  free_block(B);
  optinfo.iteration += 1;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Iteration",
      (char *) &(optinfo.iteration),sizeof(int));
  close_PSIF();
  delete [] djunk;
  delete [] disp_label;
  return 0;
}


/*** FCONST_INIT -- make sure there are _some_ force constants in PSIF_OPTKING
  1) Confirm PSIF_OPTKING has them
  2) read them from fconst.dat
  3) generate empirical force constants   ***/

void fconst_init(cartesians &carts, internals &simples, salc_set &symm) {
  int i, j, dim, count, constants_in_PSIF;
  char *buffer;
  buffer = new char[MAX_LINELENGTH];

  open_PSIF();
  constants_in_PSIF=1;
  if (psio_tocscan(PSIF_OPTKING, "Force Constants") == NULL)
    constants_in_PSIF = 0;
  close_PSIF();
  if (!constants_in_PSIF) {
    ffile_noexit(&fp_fconst, "fconst.dat",2);
    if (fp_fconst == NULL) {
      fprintf(outfile, "\nGenerating empirical Hessian.\n");
      empirical_H(simples,symm,carts);
    }
    else {
      /*** transfer force constants from fconst.dat to PSIF_OPTKING ***/
      double **temp_mat;
      dim = symm.get_num();
      temp_mat = block_matrix(dim,dim);
      for (i=0;i<dim;++i) {
        fgets(buffer,MAX_LINELENGTH,fp_fconst);
        count = 0;
        for (j=0;j<=i;++j) {
          if ( div_int(j,8) ) {
            fgets(buffer,MAX_LINELENGTH,fp_fconst);
            count = 0;
          }
          if (sscanf(buffer+count,"%lf",&(temp_mat[i][j])) != 1) {
            fprintf(outfile,"\nProblem reading force constants.\n");
            exit(2);
          }
          count += 10;
        }
      }
      fclose(fp_fconst);
      for (i=0;i<dim;++i)
        for (j=0;j<i;++j)
          temp_mat[j][i] = temp_mat[i][j];
      /*** write to PSIF_OPTKING ***/
      open_PSIF();
      psio_write_entry(PSIF_OPTKING, "Force Constants",
          (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
      close_PSIF();
      free_block(temp_mat);
    }
  }
  delete [] buffer;
}
