#include<stdio.h>
#include<math.h>
#include <ip_libv1.h>
#include <libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"init_unf_prim_atomic_grid.h"
#include"init_prun_prim_atomic_grid.h"
#include"init_unf_conc_grid.h"
#include"init_prun_conc_grid.h"
#include"free_grid_structs.h"
#include"bragg.h"
#include"SG1.h"

void grid_init(){
    int errcod,i,j;
    int num_ang_grids;
    int rpointstmp;
    int angpoints;
    int chunktmp;
    int start,end;
    int chunk_size;
    int u_atom_num;
    int Z_nuc;
    
    char *gridstring;
    
    /* Contcrete Object */
    grid_t *grid;    
    
    /* ---------------------------
       Do all of the input parsing
       and then initialize primitive
       atomic grid
       --------------------------*/
   
    grid = &(DFT_options.grid);
    
    errcod = ip_string("GRID_SPEC",&gridstring,0);
    if(errcod != IPE_OK){
	
	/* Uniform Grid */	
	
	rpointstmp = 50;
	errcod = ip_data("GRID_RPOINTS","%d",&rpointstmp,0);
	printf("rpoints = %d",rpointstmp);
	
	angpoints=302;
	errcod = ip_data("GRID_ANGPOINTS","%d",&angpoints,0);
	printf("angpoints = %d",angpoints);
	/* -------------------------
	   Construct Primitive atomic grid data 
	   ------------------------*/
	chunktmp = 4;
	errcod = ip_data("CHUNK_NUM","%d",&chunktmp,0);
	
	grid->prim_atomic_grid = init_uniform_prim_atomic_grid(rpointstmp,angpoints,chunktmp);
	init_uniform_con_grid();
       
	grid->label = "Euler-Mclaren / Lebedev Spheres";
	grid->n_rad_points = rpointstmp+2;
	for(i=0;i<Symmetry.num_unique_atoms;i++){
	    grid->atomic_grid[i].Bragg_radii = Bragg_radii[(int) Molecule.centers[Symmetry.ua2a[i]].Z_nuc];
	}
	
    }
    
    else if(!strcmp(gridstring,"SG1") || !strcmp(gridstring,"SG-1")){
       
	rpointstmp = 50;
	
	chunktmp = 4;
	errcod = ip_data("CHUNK_NUM","%d",&chunktmp,0);
	
/* Build SG1 information structures */
	grid->pruned_info.n_param_sets = 3;
	grid->pruned_info.n_diff_ang_grids = 1;
	grid->pruned_info.n_tot_ang_grids = 5;
	grid->pruned_info.a2param = (int *)init_array(Symmetry.num_unique_atoms);
	
	grid->pruned_info.param_set = (struct param_set_s *)
	    malloc(sizeof(struct param_set_s)*grid->pruned_info.n_param_sets);
	
	for(i=0;i<grid->pruned_info.n_param_sets;i++){
	    grid->pruned_info.param_set[i].n_ang_grids = 5;
	    
	    grid->pruned_info.param_set[i].alpha = 
		init_array(grid->pruned_info.param_set[i].n_ang_grids-1);
	    
	    grid->pruned_info.param_set[i].angpoints = 
		(int *) init_array(grid->pruned_info.param_set[i].n_ang_grids);
	   
	    grid->pruned_info.param_set[i].alpha[0] = SG1alpha1[i];
	    grid->pruned_info.param_set[i].alpha[1] = SG1alpha2[i];
	    grid->pruned_info.param_set[i].alpha[2] = SG1alpha3[i];
	    grid->pruned_info.param_set[i].alpha[3] = SG1alpha4[i];
	    for(j=0;j<5;j++){
		grid->pruned_info.param_set[i].angpoints[j] = SG1angular[j];
	    }
	}
	
	for(i=0;i<Symmetry.num_unique_atoms;i++){
	    Z_nuc = (int) Molecule.centers[Symmetry.ua2a[i]].Z_nuc;
	    grid->pruned_info.a2param[i] = SG1a2param[Z_nuc];
	}
	
	grid->prim_pruned_atomic_grids = init_pruned_prim_atomic_grid(rpointstmp,chunktmp,grid->pruned_info);
	init_pruned_con_grid();
	
	
	grid->label = "SG-1";
	grid->n_rad_points = 50;
	
/* Assign Bragg Radii */

	for(i=0;i<Symmetry.num_unique_atoms;i++){
	    grid->atomic_grid[i].Bragg_radii = Bragg_radii[(int) Molecule.centers[Symmetry.ua2a[i]].Z_nuc];
	}
    }
    else
	punt("No Special Grids have been implemented in this code yet");
    
    return;
}
