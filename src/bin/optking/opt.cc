/*****************************************************************************

     OPT.CC
        written by Rollin King, 1999

     main well as several auxillary functions for optking


*****************************************************************************/
extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <ctype.h>
  #include <math.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
  #include <physconst.h>
  int **get_char_table(char *symmetry);
  char **get_symm_ops(char *symmetry);
  int *get_ops_in_class(char *ptgrp);
}

#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

void intro();
void print_mat2(double **matrix, int rows, int cols, FILE *of);
double *compute_q(internals &simples, salc_set &symm);
double **compute_B(int num_atoms, internals &simples, salc_set &symm);
void delocalize(int num_atoms, internals &simples);
double **compute_G(double **B, int num_intcos, cartesians &carts);
void new_geom(cartesians &carts, internals &simples, salc_set &symm, double *dq,
	      int print_to_geom_file, int restart_geom_file, 
              char *disp_label, int disp_num, int last_disp);
void get_optinfo();
void get_syminfo(internals &simples);
void empirical_H(internals &simples, salc_set &symm, cartesians &carts);
double **compute_H(salc_set &symm,double *q, double *f_q, double **P);
double power(double x, int y);
void swap_tors(int *a, int *b, int *c, int *d);
void swap(int *a, int *b);
int *count_internals(cartesians &cart_geom, int intco_given);

int main(void)
{
  int i,j,a,b,count,dim_carts,intco_exists = 0,print_flag;
  double **B, **G, **G2, **G_inv, **H_inv, **temp_mat, **u;
  double *temp_arr2, *temp_arr, *masses;
  double *f, *f_q, *dq, *q, tval, tval2, scale, temp;
  char *disp_label; disp_label = new char[MAX_LINELENGTH];
  char buffer[MAX_LINELENGTH], *err;
  FILE *fp_fconst;

  /* These two files stay open throughout */
  ffile(&outfile, "output.dat",1) ;
  ffile(&fp_input,"input.dat",2);

  tstart(outfile);
  intro();
  get_optinfo();

  cartesians carts;
  
  dim_carts = carts.get_num_atoms()*3;
  fprintf(outfile,
    "\nGeometry and gradient from file11.dat in a.u. with masses\n");
  carts.print(2,outfile,0,disp_label,0);



  
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Generate simple internals if not given in intco.dat
    
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

  i = 0;
  fp_intco = fopen("intco.dat","r");
  if (fp_intco != NULL) {
     ip_set_uppercase(1) ;
     ip_initialize(fp_intco, outfile) ;
     if (ip_exist(":INTCO",0)) i = 1;
     ip_done();
     fclose(fp_intco);
  }

/* need to count internals coordinates in intco or compute number that will be automatically generated
   so we know how much memory to allocate (easiest way to dynamically allocate memory without rewriting lots of code) */
  number_internals = init_int_array(4);
  number_internals = count_internals(carts,i); 
    
  
  internals simples(carts,i,number_internals);
  simples.compute_internals(carts.get_num_atoms(),carts.get_coord());
  simples.compute_s(carts.get_num_atoms(),carts.get_coord() );
  fprintf(outfile,"\nSimple Internal Coordinates and Values\n");
  simples.print(outfile,1);
  fflush(outfile);

  /* Some of SYMinfo is dependent on the internal coordinates */
  get_syminfo(simples);



  
  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Generate delocalized internal coordinates if not present in intco.dat
    
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
  
  i = 1;
  fp_intco = fopen("intco.dat", "r");
  ip_set_uppercase(1) ;
  ip_initialize(fp_intco, outfile) ;
  ip_cwk_add(":INTCO");
  if (ip_exist("SYMM",0)) i = 0;
  ip_done();
  fclose(fp_intco);
  if (i) {
    if (optinfo.delocalize) {
      fprintf(outfile,"\nForming delocalized internal coordinates.\n");
      delocalize(carts.get_num_atoms(),simples);
    }
    else {
     // Print out simples to symm vector
      fp_intco = fopen("intco.dat", "r+");
      count = 0;
      for( ; ; ) {
        err = fgets(buffer, MAX_LINELENGTH, fp_intco);
        if (err == NULL) break;
        ++count;
      }

      rewind(fp_intco);
      for(j=0;j<(count-1);++j)
        err = fgets(buffer, MAX_LINELENGTH, fp_intco);
      fflush(fp_intco);
      fprintf(fp_intco,"  symm = (\n");
      for (j=0;j<simples.get_num();++j)
        fprintf(fp_intco,"    (\" \" (%d))\n",simples.index_to_id(j));
      fprintf(fp_intco,"  )\n");
      fprintf(fp_intco,")\n");
      fclose(fp_intco);
    }
  }
  fflush(outfile);




  /*^^^^^^^^^^^^^^^^^^^^^^^^^^
    Compute new geometry
    
    ^^^^^^^^^^^^^^^^^^^^^^^^*/

  if (optinfo.optimize) {
    salc_set symm("SYMM");
    symm.print();
    if (symm.get_num() == 0) {
      fprintf(outfile,"No symmetric internal coordinates to optimize.\n");
      exit(2);
    }

  // get array of values of symm internal coordinates
    dq = init_array(symm.get_num());
    q = init_array(symm.get_num());
    q = compute_q(simples,symm);

  // build G = BuB^t
    B = init_matrix(symm.get_num(),3*carts.get_num_atoms()); 
    B = compute_B(carts.get_num_atoms(),simples,symm);
    
    G = init_matrix(symm.get_num(),symm.get_num());
    G = compute_G(B,symm.get_num(),carts);

  // compute G_inv
    fprintf(outfile,"\nBuB^t ");
    G_inv = init_matrix(symm.get_num(),symm.get_num());
    G_inv = symm_matrix_invert(G,symm.get_num(),1,optinfo.redundant);
  
  // setup the masses matrix, u
    masses = init_array(3*carts.get_num_atoms());
    masses = carts.get_mass();
    u = init_matrix(dim_carts,dim_carts);
    for (i=0;i<3*carts.get_num_atoms();++i)
       u[i][i] = 1.0/masses[i];
    free(masses);

  // get forces array in cartesian coordinates, f, (in aJ/Ang)
    f = carts.forces();
   
  // compute forces in internal coordinates, f_q = G_inv B u f
    f_q = init_array(symm.get_num());
    temp_arr = init_array(symm.get_num());
    temp_arr2 = init_array(dim_carts);

    mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    mmult(B,0,&temp_arr2,1,&temp_arr,1,symm.get_num(),dim_carts,1,0);
    mmult(G_inv,0,&temp_arr,1,&f_q,1,symm.get_num(),symm.get_num(),1,0);

    free(temp_arr);
    free(temp_arr2);
    free_matrix(u,dim_carts);
   
  // Setup projection matrix P = G * G-
  // for inversion of Hessian with redundant internals
    double **P;
    P = init_matrix(symm.get_num(),symm.get_num());
    mmult(G,0,G_inv,0,P,0,symm.get_num(),symm.get_num(),symm.get_num(),0); 

  // Make Hessian if fconst.dat is absent
    fp_fconst = fopen("fconst.dat","r");
    i=0;
    if (fp_fconst == NULL)
	i=1;
    else
	fclose(fp_fconst);

    if (i) {
       fprintf(outfile, "\nGenerating empirical Hessian.\n");
       empirical_H(simples,symm,carts);
    }

  // Read in Hessian and update it if necessary from opt.aux
    H_inv = compute_H(symm,q,f_q,P);

  // Write Values and Forces of internals to opt.aux for later
    if (optinfo.bfgs) {
      ffile(&fp_opt_aux,"opt.aux",0);
      fprintf(fp_opt_aux,"Values and Forces in Intcos for Last Step.\n");
      for (i=0;i<symm.get_num();++i)
        fprintf(fp_opt_aux,"%22.16lf%22.16lf\n",q[i],f_q[i]);
      fclose(fp_opt_aux);
    }

  // Computing internal coordinate displacements H_inv f = dq, and new q
    mmult(H_inv,0,&f_q,1,&dq,1,symm.get_num(),symm.get_num(),1,0);

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
    
    free_matrix(H_inv,symm.get_num());

  // Print step summary to output.dat
    fprintf(outfile,
     "\nInternal Coordinate Update in Ang or Rad, aJ/Ang or aJ/Rad\n");
    fprintf(outfile,
     "         Value          Force          Displacement   New Value\n");
    for (i=0;i<symm.get_num();++i)
      fprintf(outfile,"%2d%15.6lf%15.6lf%15.6lf%15.6lf\n",
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
       fprintf(outfile,"\nMAX force is < %5.1e.  Optimization is complete.\n",optinfo.conv);
       fprintf(stderr,"\n  OPTKING:  optimization is complete\n");
       tstop(outfile);
       fclose(fp_input);
       fclose(outfile);
       exit(1);
    }
    free(f_q);

    strcpy(disp_label,"New Cartesian Geometry in a.u.");
    new_geom(carts,simples,symm,dq,0,0,disp_label,0,0);
    free(q); free(dq);
    free_matrix(B,symm.get_num());
    free_matrix(G_inv,symm.get_num());
  }
  else { // Doing displacements
    int num_disps = 0, disp_length = 0, restart_geom_file, line_length_count;
    double disp = 0, **displacements;

    salc_set all_salcs; // get symm and asymm salcs
    all_salcs.print();
    dq = init_array(all_salcs.get_num());

    rewind(fp_input);
    ip_initialize(fp_input, outfile) ;
    ip_set_uppercase(1) ;
    ip_cwk_add(":OPTKING") ;
    ip_count("DISPLACEMENTS",&num_disps,0);

    if (num_disps > 0) {
      displacements = init_matrix(num_disps,all_salcs.get_num());

      for (i=0;i<num_disps;++i) {
        disp_length = 0;
        ip_count("DISPLACEMENTS",&disp_length,1,i);
        if (div_int(disp_length,2) && disp_length != 0) {
          for (j=0;j<disp_length;j+=2) {
            ip_data("DISPLACEMENTS","%d",&a,2,i,j);
            ip_data("DISPLACEMENTS","%lf",&disp,2,i,j+1);
            if ((a>0) && (a<=all_salcs.get_num()))
              displacements[i][a-1] = disp; 
            else {
              punt("internal to be displaced does not exist");
            }
          }
        }
        else {
           punt("displacement vector has wrong number of elements");
          exit(1);
        }
      }
      ip_done();

      fprintf(outfile,"\nDisplaced geometries in a.u.\n");
      char *ch_temp;
      for (i=0;i<num_disps;++i)  {
        sprintf(disp_label, "Disp:");
        ch_temp = disp_label + 5 ;
        line_length_count = 5;
        for (j=0; j<all_salcs.get_num(); ++j) {
          if (fabs(displacements[i][j]) > 1.0E-8) {
            if (line_length_count < (MAX_LINELENGTH - 10)) {
              sprintf(ch_temp, " (%d %5.3lf)", j+1, displacements[i][j]);
              ch_temp += 10;
              line_length_count += 10;
            }
          }
        }
        restart_geom_file = 0;
        if (i == 0) restart_geom_file = 1;
        new_geom(carts,simples,all_salcs,displacements[i],1,
                 restart_geom_file,disp_label,i, (i==(num_disps-1)?1:0));
      }
      free_matrix(displacements,num_disps);
    }
    else {
      ip_done();
      fprintf(outfile,"\nNo DISPLACEMENTS vector found in input\n");
      exit(2);
    }
  }
  free(dq);
  tstop(outfile);
  fclose(fp_input);
  fclose(outfile);
  exit(0);
}


//Does little go into big?
int div_int(int big, int little) {
  if (little > big) return 0;
  for(;;) {
    big -= little;
    if (big == 0) return 1;
    if (big < 0) return 0;
  }
}





/*----------------------------------------------------------------------------
       PRINT_MAT2

       prints a matrix to output file

----------------------------------------------------------------------------*/

void print_mat2(double **matrix, int rows, int cols, FILE *of) {
  int i,j,col;
  for (i=0;i<rows;++i) {
    col = 0;
    for (j=0;j<cols;++j) {
      if (col == 8) {
         fprintf(outfile,"\n");
         col = 0;
      }
      fprintf(of,"%12.6f",matrix[i][j]);
      ++col;
    }
    fprintf(outfile,"\n\n");
  }
  return;
}





/*-----------------------------------------------------------------------------

	CROSS_PRODUCT

        computes cross product of two vectors
-----------------------------------------------------------------------------*/

void cross_product(double *u,double *v,double *out)
{
  out[0] = u[1]*v[2]-u[2]*v[1];
  out[1] = -1.0*(u[0]*v[2]-u[2]*v[0]);
  out[2] = u[0]*v[1]-u[1]*v[0];
  return;
}





/*----------------------------------------------------------------------------

       SCALAR_MULT

       performs scalar multiplication of a vector

----------------------------------------------------------------------------*/

void scalar_mult(double a, double *vect, int dim) {
  int i;
  for (i=0;i<dim;++i)
    vect[i] *= a;
  return;
}





/*----------------------------------------------------------------------------

       SCALAR_DIV

       performs scalar division of a vector

-----------------------------------------------------------------------------*/

void scalar_div(double a, double *vect) {
  int i;
  for (i=0;i<3;++i)
    vect[i] /= a;
  return;
}





/*----------------------------------------------------------------------------

       INTRO

       prints into

-----------------------------------------------------------------------------*/

void intro() {
fprintf(outfile,"\t_____________________________________________________\n");
fprintf(outfile,"\t|                                                   |\n");
fprintf(outfile,"\t|   OPTKING:  The Geometry Optimization Program     |\n");
fprintf(outfile,"\t|              By Rollin King, 1999                 |\n");
fprintf(outfile,"\t|___________________________________________________|\n");
}





/*----------------------------------------------------------------------------

       **SYM_MATRIX_INVERT

       inverts a matrix by diagonalization

       parameters:
             **A = matrix to be inverted
             dim = dimension of A
             print_det = print determinant if 1, nothing if 0
             redundant = zero eigenvalues allowed if 1

       returns:
             **inv_A = inverse of A

----------------------------------------------------------------------------*/

double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

  A_inv   = init_matrix(dim,dim);
  A_temp  = init_matrix(dim,dim);
  A_vects = init_matrix(dim,dim);
  A_vals  = init_array(dim);

  sq_rsp(dim,dim,A,A_vals,1,A_vects,EVAL_TOL);

  if (redundant == 0) {
     for (i=0;i<dim;++i) {
        det *= A_vals[i];
        A_inv[i][i] = 1.0/A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
     if (fabs(det) < 1E-10) {
        fprintf(outfile,"Determinant: %10.6e\n",det);
        fprintf(outfile,"Determinant is too small...aborting.\n");
        fclose(outfile);
        exit(2);
     }
  }
  else {
     for (i=0;i<dim;++i) {
        det *= A_vals[i];
        if (fabs(A_vals[i]) > REDUNDANT_EVAL_TOL)
           A_inv[i][i] = 1.0/A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

  free(A_vals);
  free_matrix(A_vects,dim);
  free_matrix(A_temp,dim);
  return A_inv;
}





/*----------------------------------------------------------------------------

	GET_OPTINFO

	reads optimization parameters from input.dat
-----------------------------------------------------------------------------*/
                                                                              
void get_optinfo() {
  int a;
  
  rewind(fp_input);
  ip_set_uppercase(1);
  ip_initialize(fp_input, outfile);
  ip_cwk_clear();
  ip_cwk_add(":OPTKING");

  /* print options */
  optinfo.print_simples = 0;
  ip_boolean("PRINT_SIMPLES", &(optinfo.print_simples),0);
  optinfo.print_params = 0;
  ip_boolean("PRINT_PARAMS", &(optinfo.print_params),0);
  optinfo.print_delocalize = 0;
  ip_boolean("PRINT_DELOCALIZE", &(optinfo.print_delocalize),0);
  optinfo.print_symmetry = 0;
  ip_boolean("PRINT_SYMMETRY", &(optinfo.print_symmetry),0);

  /* optimization parameters */
  optinfo.optimize = 1;
  if (ip_exist("DISPLACEMENTS",0)) optinfo.optimize = 0;
  optinfo.redundant = 0;
  ip_boolean("REDUNDANT", &(optinfo.redundant),0);
  optinfo.bfgs = 1;
  ip_boolean("BFGS",&(optinfo.bfgs),0);
  optinfo.mix_types = 1;
  ip_boolean("MIX_TYPES",&(optinfo.mix_types),0);
  optinfo.delocalize = 1;
  ip_boolean("DELOCALIZE",&(optinfo.delocalize),0);

  a = 5;
  ip_data("CONV","%d",&a,0);
  optinfo.conv = power(10.0, -1*a);

  a= 5;
  ip_data("EV_TOL","%d",&a,0);
  optinfo.ev_tol = power(10.0, -1*a);

  /* back-transformation parameters */
  optinfo.bt_max_iter = 500;
  ip_data("BT_MAX_ITER","%d",&(optinfo.bt_max_iter),0);
  a = 11;
  ip_data("BT_DQ_CONV","%d",&a,0);
  optinfo.bt_dq_conv  = power(10.0, -1*a);
  a = 11;
  ip_data("BT_DX_CONV","%d",&a,0);
  optinfo.bt_dx_conv  = power(10.0, -1*a);

/* Obscure limits in intco evaluation */
  optinfo.cos_tors_near_1_tol = 0.99999;
    ip_data("COS_TORS_NEAR_1_TOL","%lf",&(optinfo.cos_tors_near_1_tol),0);
  optinfo.cos_tors_near_neg1_tol = -0.99999;
    ip_data("COS_TORS_NEAR_NEG1_TOL","%lf",&(optinfo.cos_tors_near_neg1_tol),0);
  optinfo.sin_phi_denominator_tol = 0.0001;
    ip_data("SIN_PHI_DENOMINATOR_TOL","%lf",&(optinfo.sin_phi_denominator_tol),0);
    
  ip_done();
  if (optinfo.print_params) {
    fprintf(outfile,"\n+++ Optinfo Parameters +++\n");
    fprintf(outfile,"print_params:  %d\n",optinfo.print_params);
    fprintf(outfile,"print_simples: %d\n",optinfo.print_simples);
    fprintf(outfile,"print_delocalize %d\n",optinfo.print_delocalize);
    fprintf(outfile,"print_symmetry %d\n",optinfo.print_symmetry);
    fprintf(outfile,"optimize:      %d\n",optinfo.optimize);
    fprintf(outfile,"redundant:     %d\n",optinfo.redundant);
    fprintf(outfile,"bfgs:          %d\n",optinfo.bfgs);
    fprintf(outfile,"mix_types:     %d\n",optinfo.mix_types);
    fprintf(outfile,"conv:          %.1e\n",optinfo.conv);
    fprintf(outfile,"bt_max_iter:   %d\n",optinfo.bt_max_iter);
    fprintf(outfile,"bt_max_dq_conv:    %.1e\n",optinfo.bt_dq_conv);
    fprintf(outfile,"bt_max_dx_conv:    %.1e\n",optinfo.bt_dx_conv);
    fprintf(outfile,"cos_tors_near_1_tol:     %10.6e\n",
        optinfo.cos_tors_near_1_tol);
    fprintf(outfile,"cos_tors_near_neg1_tol: %10.6e\n",
        optinfo.cos_tors_near_neg1_tol);
    fprintf(outfile,"sin_phi_denominator_tol: %10.6e\n",
        optinfo.sin_phi_denominator_tol);
    fflush(outfile);
  }
}





/*------------------------------------------------------------------------------

       GET_SYMINFO

       gets symmetry info
-----------------------------------------------------------------------------*/ 

void get_syminfo(internals &simples) {
  int a, b, c, d, aa, bb, cc, dd, i, j, sign;
  int id, intco_type, sub_index, ops, num_atoms;

  rewind(fp_input);
  ip_set_uppercase(1);
  ip_initialize(fp_input,outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  file30_init();

  if (ip_exist("SYMMETRY",0))
    ip_string("SYMMETRY",&syminfo.symmetry,0);
  else {
    syminfo.symmetry = file30_rd_sym_label();
  }
  num_irreps = syminfo.num_irreps = file30_rd_nirreps();
  syminfo.ict = file30_rd_ict();
  syminfo.irrep_lbls = file30_rd_irr_labs();
  num_atoms = file30_rd_natom();
  file30_close();
  ip_done();

  j = strlen(syminfo.symmetry);
  strncpy(ptgrp,syminfo.symmetry,j);
  for ( ;j<3;++j)
    ptgrp[j] = ' ';
  ptgrp[3] = '\0';
  for (i=0;i<3;++i) 
      ptgrp[i] = toupper(ptgrp[i]);

  syminfo.clean_irrep_lbls = new char*[num_irreps];
  for (i=0;i<num_irreps;++i) {
    syminfo.clean_irrep_lbls[i] = new char[4];
    strcpy(syminfo.clean_irrep_lbls[i], syminfo.irrep_lbls[i]);
    for (j=0;j<3;++j) {
      if (syminfo.clean_irrep_lbls[i][j] == '\"')
         syminfo.clean_irrep_lbls[i][j] = 'P' ;
      if (syminfo.clean_irrep_lbls[i][j] == '\'')
         syminfo.clean_irrep_lbls[i][j] = 'p' ;
      }
    }
  syminfo.ct = get_char_table(ptgrp);
  syminfo.op_lbls = get_symm_ops(ptgrp);
  syminfo.ict_ops = init_int_matrix(simples.get_num(),num_irreps);
  syminfo.ict_ops_sign = init_int_matrix(simples.get_num(),num_irreps);
  ops_in_class = init_int_array(num_irreps);
  ops_in_class = get_ops_in_class(ptgrp);

  for (i=0;i<simples.get_num();++i)
    for (j=0;j<num_irreps;++j)
      syminfo.ict_ops_sign[i][j] = 1;

 // Generate simple internal coordinate transformation matrix
  for (i=0;i<simples.get_num();++i) {
     id = simples.index_to_id(i);
     simples.locate_id(id,&intco_type,&sub_index);
     if (intco_type == STRE_TYPE) {
        a = simples.stre.get_A(sub_index);
        b = simples.stre.get_B(sub_index);
        for (ops=0;ops < num_irreps;++ops) {
          aa = syminfo.ict[ops][a]-1;
          bb = syminfo.ict[ops][b]-1;
          swap(&aa,&bb);
          syminfo.ict_ops[i][ops] = simples.stre.get_id_from_atoms(aa,bb);
        }
     }
     if (intco_type == BEND_TYPE) {
        a = simples.bend.get_A(sub_index);
        b = simples.bend.get_B(sub_index);
        c = simples.bend.get_C(sub_index);
        for (ops=0;ops < num_irreps;++ops) {
          aa = syminfo.ict[ops][a]-1;
          bb = syminfo.ict[ops][b]-1;
          cc = syminfo.ict[ops][c]-1;
          swap(&aa,&cc);
          syminfo.ict_ops[i][ops] = simples.bend.get_id_from_atoms(aa,bb,cc);
        }
     }
     if (intco_type == TORS_TYPE) {
        a = simples.tors.get_A(sub_index);
        b = simples.tors.get_B(sub_index);
        c = simples.tors.get_C(sub_index);
        d = simples.tors.get_D(sub_index);
        for (ops=0;ops < num_irreps;++ops) {
          aa = syminfo.ict[ops][a]-1;
          bb = syminfo.ict[ops][b]-1;
          cc = syminfo.ict[ops][c]-1;
          dd = syminfo.ict[ops][d]-1;
          swap_tors(&aa, &bb, &cc, &dd);
          syminfo.ict_ops[i][ops] = simples.tors.get_id_from_atoms(aa,bb,cc,dd);
          if ( ('S' == syminfo.op_lbls[ops][0]) ||
               ('I' == syminfo.op_lbls[ops][0]) )
               syminfo.ict_ops_sign[i][ops] = -1;
        }
     }
     if (intco_type == OUT_TYPE) {
        a = simples.out.get_A(sub_index);
        b = simples.out.get_B(sub_index);
        c = simples.out.get_C(sub_index);
        d = simples.out.get_D(sub_index);
        for (ops=0;ops < num_irreps;++ops) {
          aa = syminfo.ict[ops][a]-1;
          bb = syminfo.ict[ops][b]-1;
          cc = syminfo.ict[ops][c]-1;
          dd = syminfo.ict[ops][d]-1;
          syminfo.ict_ops[i][ops] = simples.out.get_id_from_atoms(aa,bb,cc,dd,&sign);
          if ( ('S' == syminfo.op_lbls[ops][0]) ||
               ('I' == syminfo.op_lbls[ops][0]) )
               sign *= -1;
          syminfo.ict_ops_sign[i][ops] = sign;
        }
     }
  }

  if (optinfo.print_symmetry) {
    fprintf(outfile,"\n+++ Symmetry Information +++\n");
    fprintf(outfile,"The ICT table from file30:\n");
    for(i=0;i<num_irreps;++i) {
       for(j=0;j<num_atoms;++j)
          fprintf(outfile,"%3d",syminfo.ict[i][j]);
       fprintf(outfile,"\n");
    }

    fprintf(outfile,"Clean irrep labels:");
    for (j=0;j<num_irreps;++j)
      fprintf(outfile,"%5s ", syminfo.clean_irrep_lbls[j]);
    fprintf(outfile,"\n");

    fprintf(outfile,"\nCharacter table from char_table.c and symmetry: %s\n", ptgrp);
    fprintf(outfile,"      ");
    for (i=0;i<num_irreps;++i) fprintf(outfile,"%5s",syminfo.op_lbls[i]);
    fprintf(outfile,"\n");

    for (i=0;i<num_irreps;++i) {
      for (j=0;j<num_irreps;++j) {
        if (j == 0) fprintf(outfile,"%5s ", syminfo.irrep_lbls[i]);
        fprintf(outfile,"%5d",syminfo.ct[i][j]);
      }
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"Internal coordinate transformation matrix.\n");
    for (i=0;i<num_irreps;++i) fprintf(outfile,"%5s",syminfo.op_lbls[i]);
    fprintf(outfile,"\n");
    for (id=0;id<simples.get_num();++id) {
      for (ops=0;ops < num_irreps;++ops)
         fprintf(outfile,"%5d",syminfo.ict_ops[id][ops]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,
      "Internal transformation sign matrix to fix torsions and out-of-planes.\n");
    for (i=0;i<num_irreps;++i) fprintf(outfile,"%5s",syminfo.op_lbls[i]);
    fprintf(outfile,"\n");
    for (id=0;id<simples.get_num();++id) {
      for (ops=0;ops < num_irreps;++ops)
         fprintf(outfile,"%5d",syminfo.ict_ops_sign[id][ops]);
      fprintf(outfile,"\n");
    }
  }
  return;
}





/*------------------------------------------------------------------------------

       POWER

       raises number to a power
------------------------------------------------------------------------------*/

double power(double x, int y) {
  double tval = 1.0;
  int invert = 0;

  if (y < 0) {
     invert = 1;
     y = abs(y);
  }
  for( ; y>0 ; --y)
    tval *= x;
  if (invert) tval = 1.0/tval;
  return tval;
}





/*------------------------------------------------------------------------------

       SWAP_TORS

       swaps a and d, b and c ... something about torsional angles??
------------------------------------------------------------------------------*/

void swap_tors(int *a, int *b, int *c, int *d) {
  int p;
  if (*a > *d) {
     p = *a;
    *a = *d;
    *d = p;
     p = *b;
    *b = *c;
    *c = p;
  }
}






/*------------------------------------------------------------------------------

       SWAP

       swaps a <-> b
-----------------------------------------------------------------------------*/

void swap(int *a, int *b) {
  int c;
  if (*a > *b) {
    c = *a;
    *a = *b;
    *b = c;
  }
  return;
}




/*------------------------------------------------------------------------------

       COUNT_INTERNALS

       purpose: if intco contains internal coordinates count them
                else count number of internals which will be automatically
		generated

       parameters:
            &cart_geom -- address of object of type cartesian
	    intco_given -- 1 = internals given in intco
	                   0 = no internals in intco

       returns:
            count_array -- array containing the numbers of each type of internal

   *this is basically ripped off of the constructor for the internals class
	                                                   J. Kenny 6/30/00
------------------------------------------------------------------------------*/

int *count_internals(cartesians &cart_geom, int intco_given) {

  int i,j,k,l,a,b,
      num, num_atoms, count,      /* counter variables */
      Zmax,Zmin,
      *ioff,                      /* ioff array for indexing */
      *count_array,               /* array to hold number of each internal type */
      **bonds;                    /* 1 if atoms bonded, 0 otherwise */

  double *coord,                  /* holds coordinates from &cart_geom */
         *distance;               /* holds computed distances */

  count_array = init_int_array(4);
      
    /*#########################################
      if intco_given=1 count internals in intco
    #########################################*/
  if (intco_given) {

      /* initialize intco.dat */
      ip_done();
      fp_intco = fopen("intco.dat","r");
      ip_set_uppercase(1);
      ip_initialize(fp_intco, outfile);
      ip_cwk_add(":INTCO");

      /* count stretches */
      num=0;
      if (ip_exist("STRE",0)) {
	  ip_count("STRE",&num,0);
	}
      count_array[0]=num;

      /* count bends */
      num=0;
      if (ip_exist("BEND",0)) {
	  ip_count("BEND",&num,0);
	}
      count_array[1]=num;

      /* count torsions */
      num=0;
      if (ip_exist("TORS",0)) {
	  ip_count("TORS",&num,0);
	}
      count_array[2]=num;

      /* count out of plane bends */
      num=0;
      if (ip_exist("OUT",0)) {
	  ip_count("OUT",&num,0);
	}
      count_array[3]=num;

      fclose(fp_intco);
      ip_done();
    }

  
  /*####################################################################
    if no internals in intco.dat count the number that will be generated
   ###################################################################*/
  else {

      num_atoms = cart_geom.get_num_atoms();
      coord = init_array(3 * num_atoms);
      coord = cart_geom.get_coord();
   

      /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	figure out which atoms are bonded (I'm not exactly sure how this works)
	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      ioff = (int *) malloc (32641 * sizeof(int));
      ioff[0]=0;
      for (i=1;i<32641;++i) {
	  ioff[i] = ioff[i-1] + i;
	}
		  
      /* compute distances */
      distance = init_array( ((num_atoms+1)*num_atoms)/2 );
      count = -1;
      for(i=0;i<num_atoms;++i) {
	  for(j=0;j<=i;++j) {
	      distance[++count] = sqrt(SQR(coord[3*i+0]-coord[3*j+0]) +
				       SQR(coord[3*i+1]-coord[3*j+1]) +
				       SQR(coord[3*i+2]-coord[3*j+2]));
	    }
	}

      /* determine bonds */
      bonds = init_int_matrix(num_atoms,num_atoms);
      for (i=0;i<num_atoms;++i) {
	  for(j=0;j<=i;++j) {
	      Zmax = MAX((int)cart_geom.get_atomic_num(i),(int)cart_geom.get_atomic_num(j));
	      Zmin = MIN((int)cart_geom.get_atomic_num(i),(int)cart_geom.get_atomic_num(j));
	      a = ioff[Zmax-1] + (Zmin-1);
	      if (bondl[a] != 0.0) {
		  if (distance[ioff[i]+j] < (1.2 * bondl[a])) {
		      bonds[i][j] = 1;
		      bonds[j][i] = 1;
		    }
		}
	      else {
		  fprintf(outfile,"WARNING! Optking does not know what bond lengths");
		  fprintf(outfile,"to expect for all the atoms.\n");
		  fprintf(outfile,"You may have to specify connectivity in input.");
		}
	    }
	}


      /* check input for user specified bonds or nobonds */
      rewind(fp_input);
      ip_set_uppercase(1);
      ip_initialize(fp_input,outfile);
      ip_cwk_add(":OPTKING");

      if (ip_exist("BONDS",0)) {
	  ip_count("BONDS",&num,0);
	  for(i=0;i<num;++i) {
	      ip_data("BONDS","%d",&a,2,i,0);
	      ip_data("BONDS","%d",&b,2,i,0);
	      a -= 1;  b -= 1;
	      bonds[a][b] = 1;
	      bonds[b][a] = 1;
	    }
	}

       if (ip_exist("NOBONDS",0)) {
	  ip_count("NOBONDS",&num,0);
	  for(i=0;i<num;++i) {
	      ip_data("NOBONDS","%d",&a,2,i,0);
	      ip_data("NOBONDS","%d",&b,2,i,0);
	      a -= 1;  b -= 1;
	      bonds[a][b] = 0;
	      bonds[b][a] = 0;
	    }
	 }

       ip_done();

	/* count number of bonds */
	num=0;
	for(i=0;i<num_atoms;++i) {
	    for(j=i+1;j<num_atoms;++j) {
		if(bonds[i][j] == 1) {
		    ++num;
		  }
	      }
	  }
	count_array[0]=num;


	/*^^^^^^^^^^^^^^^^^^^^^^^
	  count number of bends
	  ^^^^^^^^^^^^^^^^^^^^^*/
        num=0;
	for(i=0;i<num_atoms;++i) {
	    for(j=0;j<num_atoms;++j) {
		if(i!=j)
		for(k=i+1;k<num_atoms;++k) {
		    if(j!=k)
		    if (bonds[i][j] && bonds[j][k]) {
			++num;
		      }
		  }
	      }
	  }
         count_array[1]=num;


	 /*^^^^^^^^^^^^^^^^^^^^^^^^^
	   count number of torsions
	   ^^^^^^^^^^^^^^^^^^^^^^^*/
         num=0;
	 for(i=0;i<num_atoms;++i) {
	     for(j=0;j<num_atoms;++j) {
		 if(i!=j)
		 for(k=0;k<num_atoms;++k) {
		     if(i != k && j != k) {
			 for(l=i+1;l<num_atoms;++l) {
			     if( (l != j && l != k) && bonds[i][j] && bonds[j][k] && bonds[k][l]) {
				 ++num;
			       }
			   }
		       }
		   }
	       }
	   }
	 count_array[2]=num;
	 count_array[3]=0;

         free(ioff);
         free(distance);
         free(coord);
         free_int_matrix(bonds,num_atoms);
    }

   return count_array;

}	 
		
	       


/*-------------------
  PUNT()

  prints errors
  and exits
  -----------------*/
extern "C" {
void punt( char *message ) {
  fprintf(outfile,"\nerror: %s\n", message);
  fprintf(outfile,"         *** stopping execution ***\n");
  fprintf(stderr,"\n OPTKING error: %s\n", message);
  fprintf(stderr,"                 *** stopping execution ***\n");
  fclose(outfile);
  exit(1);
}
}





