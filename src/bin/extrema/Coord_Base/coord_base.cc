/*###########################################################################*/
/*! \file coord_base.cc
  \brief Actual constructor and member functions for coordinate base class. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#define EXTERN
#include "extrema.h"



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::construct_coord_base(int natm, int nent, int nopt) 
  \brief Actual constructor.
  
  This is the actual constructor.  The parameters may not be known when 
  a derived class is initialized and this function may be called 
  later in the derived class's constructor.  This obviously undermines
  the idea of automatic construction and developers must be careful to 
  call this function at the appropriate point.
  \param natm number of atoms
  \param nent number of entries (atoms + dummy atoms)
  \param nopt number of optimized coordinates */
/*---------------------------------------------------------------------------*/

void coord_base :: construct_coord_base(int natm, int nent, int nopt) {

    num_atoms = natm;
    num_entries = nent;
    num_coords =  nopt;

    parse_input();
    
    carts = init_array(3*num_entries);
    c_grads = init_array(3*num_entries);
    coords = init_array(num_coords);
    grads = init_array(num_coords);
    Hi = init_matrix(num_coords,num_coords);
    u = init_matrix(3*num_entries,3*num_entries);
    masses = init_array(num_entries);
    coords_old = init_array(num_coords);
    grads_old = init_array(num_coords);  
    Hi_old = init_matrix(num_coords,num_coords);
    coord_write = init_array(num_coords);
   
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::update_Hi()
  \brief Updates inverse hessian.

  Interface for <b>math_tools</b> update functions. */
/*---------------------------------------------------------------------------*/

void coord_base :: update_Hi() {

    int i;
    double *coord_dif, *grad_dif;
    coord_dif = init_array(num_coords);
    grad_dif = init_array(num_coords);
    for(i=0;i<num_coords;++i) {
	coord_dif[i] = coords[i] - coords_old[i];
	grad_dif[i] = grads[i] - grads_old[i];
    }

    if(!strcmp(update,"MS")) {
	fprintf(outfile,"\n  Performing ms update of inverse hessian\n");
	Hi = math_tools :: update_ms(num_coords, coord_dif, grad_dif, Hi_old);
    }
    else { 
	fprintf(outfile,"\n  Performing bfgs update of inverse hessian\n");
	Hi = math_tools ::update_bfgs(num_coords, coord_dif, grad_dif, Hi_old);
    }

    free(coord_dif);
    free(grad_dif);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::grad_test()
  \brief Tests for convergence of cartesian gradients. */
/*---------------------------------------------------------------------------*/

void coord_base :: grad_test() {
    
    int i, conv=1;
    double sum=0.0;

    for(i=0;i<3*num_entries;++i) 
	sum += fabs(c_grads[i]);
    
    sum /= (3*num_atoms);

    fprintf(outfile,"\n  RMS gradient = %.10lf\n", sum);

    for(i=0;i<(3*num_entries);++i) 
	if(fabs(c_grads[i]) > (1.0/pow(10.0,(double)grad_max)))
	    conv = 0;
    if(conv) {
	fprintf(outfile,"\n  All gradients below convergence criteria of");
	fprintf(outfile," 10^-%d",grad_max);
	fprintf(outfile,"\n  Optimization completed\n");
	fprintf(stdout,"\n Optimization completed");  
	 tstop(outfile);
	 fclose(infile);
	 fclose(outfile);
	 exit(1);
    }
    
    return;
}



/*--------------------------------------------------------------------------*/
/*! \fn coord_base::H_test()
  \brief Computes hessian and its eigenvalues. */
/*-------------------------------------------------------------------------*/

void coord_base :: H_test() {

    int i;
    double *evals, **evecs, **H;

    evals = init_array(num_coords);
    evecs = init_matrix(num_coords,num_coords);

    H = symm_matrix_invert(Hi,num_coords,0,0);

    sq_rsp( num_coords, num_coords, H, evals, 1, evecs, 1.0e-14); 
    
    if(print_lvl > NORMAL_PRINT) {
	fprintf(outfile,"\n  H matrix (a.u.):\n");
	print_mat(H, num_coords, num_coords, outfile);

	fprintf(outfile,"\n  H eigenvalues:");
	for(i=0;i<num_coords;++i)
	    fprintf(outfile,"\n  %lf",evals[i]);
	fprintf(outfile,"\n");
    }

    free(evals);
    free_matrix(evecs, num_coords);

    return;
}
