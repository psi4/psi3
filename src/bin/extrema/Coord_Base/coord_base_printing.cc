/*###########################################################################*/
/*! \file coord_base_printing.cc
  \brief  Printing member functions for coordinate base class. */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"



/*---------------------------------------------------------------------------*/
/*! \fn coord_base::print_carts(double)
  \brief Prints cartesians.

  \param conv conversion factor;
  either <b>1.0</b> for bohr or <b>_bohr2angstroms</b> for angstroms */
/*---------------------------------------------------------------------------*/

void coord_base :: print_carts(double conv) {

    int i, j;
    double **temp;
    
    temp = init_matrix(num_entries,3);
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    temp[i][j] = carts[3*i+j]*conv;
    if(conv==1.0)
	fprintf(outfile,"\n  Cartesian Coordinates (bohr):\n");
    else
	fprintf(outfile,"\n  Cartesian Coordinates (angstroms):\n");
     fprintf(outfile,"                       x              y         ");
    fprintf(outfile,"       z\n");
    fprintf(outfile,"                --------------- --------------- ");
    fprintf(outfile,"---------------\n");
    for(i=0;i<num_entries;++i)
	fprintf(outfile,"  %-10s  %15.10lf %15.10lf %15.10lf\n",
		felement[i], temp[i][0], temp[i][1], temp[i][2]);
    free_matrix(temp,num_entries);

    return;
}



/*--------------------------------------------------------------------------*/
/*! \fn coord_base::print_c_grads()
  \brief Prints cartesian gradients (Hartree/Bohr). */
/*--------------------------------------------------------------------------*/

void coord_base :: print_c_grads() {
	
    int i, j;
    double **temp;
    temp = init_matrix(num_entries,3);
	
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    temp[i][j] = c_grads[3*i+j];
    fprintf(outfile,"\n  Cartesian Gradients (a.u):\n");
    fprintf(outfile,"                       x              y         ");
    fprintf(outfile,"       z\n");
    fprintf(outfile,"                --------------- --------------- ");
    fprintf(outfile,"---------------\n");
    for(i=0;i<num_entries;++i)
	fprintf(outfile,"  %-10s  %15.10lf %15.10lf %15.10lf\n",
		felement[i], temp[i][0], temp[i][1], temp[i][2]);
    free_matrix(temp,num_entries);
    
    return;
}


/*---------------------------------------------------------------------------*/
/*! \fn coord_base::print_Hi()
  \brief prints inverse of hessian matrix. */
/*---------------------------------------------------------------------------*/

void coord_base :: print_Hi() {

    fprintf(outfile,"\n  Inverse hessian matrix(a.u.):\n");
    print_mat(Hi,num_coords,num_coords,outfile);

    return;
}



