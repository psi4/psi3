#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

inline double fast_pow(double x, int l);

void compute_grid_mos()
{
  int i,j,k,l,ig,jg,ibf,jbf,ib,jb,jlim,kk,ll;
  int iatom, jatom, iang, jang, i_1stbf, j_1stbf, nbfi, nbfj;
  int igmin, igmax;
  int ixm, iym, izm, iind;
  int ix, iy, iz;
  double ax, ay, az, xa, ya, za, ra2;
  double ai, norm_pf, ang_pf, exponent;
  double x,y,z;
  double r,r2,tmp;
  double *bf_values, *values;
  double *evec;


  /* Initialize some intermediates */
  grid3d_pts = init_box(nix+1,niy+1,niz+1);
  bf_values = init_array(nbfso);
  evec = init_array(nbfso);

  for(i=0;i<nbfso;i++)
      evec[i] = scf_evec_ao[i][mo_to_plot];

  for(ix=0;ix<=nix;ix++) {
      for(iy=0;iy<=niy;iy++) {
	x = grid_origin[0] + grid_step_x[0]*ix + grid_step_y[0]*iy;
	y = grid_origin[1] + grid_step_x[1]*ix + grid_step_y[1]*iy;
	z = grid_origin[2] + grid_step_x[2]*ix + grid_step_y[2]*iy;

	for(iz=0;iz<=niz;iz++) {
	  /* Shell loop */
	  for(i=0;i<nshell;i++) {
	      iatom = snuc[i] - 1;
	      ax = geom[iatom][0];
	      ay = geom[iatom][1];
	      az = geom[iatom][2];
	      iang = stype[i]-1;
	      izm = 1;
	      iym = iang+1;
	      ixm = iym*iym;
	      i_1stbf = sloc[i] - 1;
	      nbfi = (iang+2)*(iang+1)/2;
	      igmin = sprim[i] - 1;
	      igmax = igmin + snumg[i] - 1;
    
	      xa = x - ax;
	      ya = y - ay;
	      za = z - az;
	      ra2 = xa*xa + ya*ya + za*za;

	      /*--- Compute exponentional factor - loop over primitives ---*/
	      exponent = 0.0;
	      for(ig=igmin;ig<=igmax;ig++) {
		  ai = exps[ig];
		  exponent += contr[ig]*exp(-ai*ra2);
	      }

	      /*--- Loop over basis functions in the shell ---*/
	      values = &(bf_values[i_1stbf]);
	      for(ibf=0;ibf<nbfi;ibf++,values++) {
		  norm_pf = norm_bf[iang][ibf];
		  ang_pf = fast_pow(xa,xpow_bf[iang][ibf])*
			   fast_pow(ya,ypow_bf[iang][ibf])*
			   fast_pow(za,zpow_bf[iang][ibf]);
		  
		  *values = norm_pf * ang_pf * exponent;
	      }
	  } /*--- end of shell loop ---*/
	  tmp = 0;
	  for(i=0;i<nbfso;i++)
	      tmp += bf_values[i]*evec[i];
	  grid3d_pts[ix][iy][iz] = tmp;
	  
          x += grid_step_z[0];
	  y += grid_step_z[1];
	  z += grid_step_z[2];
	}
      }
  }
}

inline double fast_pow(double x, int l)
{
  switch (l) {
  case 0:
      return 1.0;
  case 1:
      return x;
  case 2:
      return x*x;
  case 3:
      return x*x*x;
  case 4:
      return (x*x)*(x*x);
  default:
      return pow(x,l);
  }
}

