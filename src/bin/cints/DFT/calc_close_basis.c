#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"bas_comp_functions.h"
#define TOL 1E-10

void calc_close_basis(int atom_num, int chunk_num){
    int i,j,k,l;
    int chunk_center;
    int shell_center;
    int am2shell;
    int max_am;
    int num_ao;
    int num_shells;
    
    double chunk_rad_in;
    double chunk_rad_out;
    double chunk_rad;
    double bragg;
    double xd,yd,zd,rr,r;
    double bastmp;
    
    struct coordinates point_geom;
    struct coordinates atom_point_geom;
    struct coordinates shell_geom;
    
    struct atomic_grid_s *atom_grid;
    struct leb_chunk_s *chunk;
    
    num_shells = BasisSet.num_shells;
    max_am = BasisSet.max_am;

    atom_grid = &(DFT_options.grid.atomic_grid[atom_num]);
    bragg = atom_grid->Bragg_radii;
    atom_point_geom = atom_grid->atom_center;
    chunk_center = atom_grid->atom_num;

    chunk = &(atom_grid->leb_chunk[chunk_num]);   
    chunk_rad_in = chunk->spheres[0].r*bragg;
    chunk_rad_out = chunk->spheres[chunk->size-1].r*bragg;

    
    
    for(i=0;i<num_shells;i++){
	
	am2shell = BasisSet.am2shell[i];
	shell_center = BasisSet.shells[am2shell].center - 1;
	shell_geom = Molecule.centers[shell_center];
	
	if(shell_center == chunk_center){
	    rr=chunk_rad_in*chunk_rad_in;
	}    
	else{
	    chunk_rad=chunk_rad_out;
	/* ---------------------------------
	   Determine where to calculate the 
	   basis function
	   --------------------------------*/
    
	    point_geom.x = chunk_rad*(shell_geom.x-atom_point_geom.x)
		+atom_point_geom.x;
	    point_geom.y = chunk_rad*(shell_geom.y-atom_point_geom.y)
		+atom_point_geom.y;
	    point_geom.z = chunk_rad*(shell_geom.z-atom_point_geom.z)
		+atom_point_geom.z;
	    
	    xd = point_geom.x-shell_geom.x;
	    yd = point_geom.y-shell_geom.y;
	    zd = point_geom.z-shell_geom.z;
	    
	    rr = xd*xd+yd*yd+zd*zd;
	    r = sqrt(rr);
	    
	    
	}
	
	bastmp = calc_radial_bas(am2shell,rr,r);
	
	fprintf(outfile,"\nchunk_num = %d shell_num = %d rr = %10.10lf bastmp = %10.10lf",
		chunk_num,i,rr,bastmp);
    }
}
	
    
    
