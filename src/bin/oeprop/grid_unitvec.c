#define EXTERN
#include "includes.h"
#include "oeprop.h"
#include "oeprop.gbl"


void grid_unitvec()
{
  int i;
  double sum1, sum2, dot;
  double step_x, step_y;

  /* Normalizing the grid_unit's */

  dot_arr(grid_unit_x,grid_unit_x,3,&sum1);
  dot_arr(grid_unit_y,grid_unit_y,3,&sum2);
  sum1 = sqrt(sum1); sum2 = sqrt(sum2);
  for(i=0;i<3;i++) {
    grid_unit_x[i] /= sum1;
    grid_unit_y[i] /= sum2;
  }
  

  /* Checking if vectors are parallel */

  dot_arr(grid_unit_x,grid_unit_y,3,&dot);
  if (1.0 - fabs(dot) < ADOTB_ORTHOGONAL) { /* Vectors are parallel - aborting */
    fprintf(outfile,"Vectors GRID_UNIT_X and GRID_UNIT_Y are parallel. Aborting.\n\n");
    exit(2);
  }
  else if (fabs(dot) > ADOTB_ORTHOGONAL) {  /* Vectors are not orthogonal - orthonormalizing */
         for(i=0;i<3;i++)
	   grid_unit_y[i] -= dot*grid_unit_x[i];
         dot_arr(grid_unit_y,grid_unit_y,3,&sum1);
	 sum1 = sqrt(sum1);
	 for(i=0;i<3;i++)
	   grid_unit_y[i] /= sum1;
       }

  
  /* Moving the origin into the lower left corner of the grid rectangle */

  for(i=0;i<3;i++)
    grid_origin[i] += grid_xy0[0]*grid_unit_x[i] + grid_xy0[1]*grid_unit_y[i];


  /* Computing grid size and unit cell vectors along x and y */

  step_x = (grid_xy1[0] - grid_xy0[0])/nix;
  step_y = (grid_xy1[1] - grid_xy0[1])/niy;
  for(i=0;i<3;i++) {
    grid_step_x[i] = grid_unit_x[i]*step_x;
    grid_step_y[i] = grid_unit_y[i]*step_y;
  }
}


