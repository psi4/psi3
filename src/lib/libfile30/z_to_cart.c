#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <physconst.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_z_to_cart(z_entry *z_geom, int num_atoms)
**  : return cartesian coordinates for given z-matrix
**
**  arguments:
**    z_entry *z_geom -- array num_atoms long containing z-matrix info
**    int num_atoms   -- number of atoms
** 
**  returns: 
**    double **cart_geom -- num_atoms x 3 matrix containing cartesian coords
**
**  NOTE: this code assumes valid z-matrix is passed to it ... no error checking
*/

double *cross_pdt(double *u, double *v);
double *unit(double *u, double *v);
double dot_pdt(double *u, double *v);

double **file30_z_to_cart( struct z_entry *z_geom, int num_atoms ) {

  int linear = 1, i, j;
  
  double linear_cutoff = 1.0E-6,
         **geom, *eAB, *ez, *ey, *ex,
	 cos_ABC, sin_ABC, cos_BCD, sin_BCD, cos_ABCD, sin_ABCD;

  geom = init_matrix(num_atoms,3);

  for(i=0;i<num_atoms;++i) {

      if(i==0) {
	  geom[0][0] = geom[0][1] = geom[0][2] = 0.0;
	}

      else if(i==1) {
	  geom[1][0] = geom[1][1] = 0.0;
	  geom[1][2] = z_geom[i].bond_val;
	}

      else if(i==2) {
	  if(z_geom[2].bond_atom == 2) {
	      geom[2][0] = geom[1][0] + z_geom[2].bond_val * sin(z_geom[2].angle_val * _pi/180.0);
	      geom[2][1] = geom[1][1];
	      geom[2][2] = geom[1][2] - z_geom[2].bond_val * cos(z_geom[2].angle_val * _pi/180.0);
	    }
	  if(z_geom[2].bond_atom == 1) {
	      geom[2][0] = geom[0][0] + z_geom[2].bond_val * sin(z_geom[2].angle_val * _pi/180.0);
	      geom[2][1] = geom[0][1];
	      geom[2][2] = geom[0][2] + z_geom[2].bond_val * cos(z_geom[2].angle_val * _pi/180.0);
	    }
	}

      else {

	  /*this section handles linear fragments in the beginning of the z-matrix*/
	  if( ((geom[i-1][0] - linear_cutoff) < 0.0) && linear ) {

	      if( (geom[z_geom[i].bond_atom][2] - geom[z_geom[i].angle_atom][2]) > 0.0 ) {
		  geom[i][0] = z_geom[i].bond_val * sin(z_geom[i].angle_val * _pi/180.0);
		  geom[i][1] = 0.0;
		  geom[i][2] = geom[z_geom[i].bond_atom-1][2] - z_geom[i].bond_val * cos(z_geom[i].angle_val * _pi/180.0);
		}

	      else {
		  geom[i][0] = z_geom[i].bond_val * sin(z_geom[i].angle_val * _pi/180.0);
		  geom[i][1] = 0.0;
		  geom[i][2] = geom[z_geom[i].bond_atom-1][2] + z_geom[i].bond_val * cos(z_geom[i].angle_val * _pi/180.0);
		}
	    }

	  /*this handles nonlinear fragments*/
	  else {

	      linear = 0;

	      eAB = unit(geom[z_geom[i].tors_atom-1],geom[z_geom[i].angle_atom-1]);
	      ez  = unit(geom[z_geom[i].angle_atom-1],geom[z_geom[i].bond_atom-1]);
              ey = cross_pdt(eAB,ez);
	      cos_ABC = -dot_pdt(ez,eAB);
              sin_ABC = sqrt(1 - (cos_ABC * cos_ABC) );
              cos_BCD = cos( z_geom[i].angle_val * _pi/180.0);
              sin_BCD = sin( z_geom[i].angle_val * _pi/180.0);
              cos_ABCD = cos( z_geom[i].tors_val * _pi/180.0);
              sin_ABCD = sin( z_geom[i].tors_val * _pi/180.0);
	      for(j=0;j<3;++j) ey[j] /= sin_ABC;
              ex = cross_pdt(ey,ez);

	      for(j=0;j<3;++j) {
		  geom[i][j] = geom[z_geom[i].bond_atom-1][j] + z_geom[i].bond_val *
				    ( -ez[j]*cos_BCD + ex[j]*sin_BCD*cos_ABCD +ey[j]*sin_BCD*sin_ABCD );
		}

	      free(eAB);
              free(ez);
              free(ey);
              free(ex);
	    } 
	}
    }

  return geom;
}
	      
		   
double *cross_pdt(double *u,double *v) {

  double *out;
  out = init_array(3);
  
  out[0] = u[1]*v[2]-u[2]*v[1];
  out[1] = -1.0*(u[0]*v[2]-u[2]*v[0]);
  out[2] = u[0]*v[1]-u[1]*v[0];

  return out;
}

double *unit(double *u, double *v) {

  int i;
  double *out, norm=0;
  out = init_array(3);

  for(i=0;i<3;++i) {
      out[i] = v[i] - u[i];
      norm += out[i] *out[i];
    }

  norm = sqrt(norm);

  for(i=0;i<3;++i) out[i] /= norm;

  return out;
}


double dot_pdt(double *u, double *v) {

  int i;
  double sum = 0;

  for(i=0;i<3;i++) {
    sum += u[i]*v[i];
  }

  return sum;
}
	  
	  
