#include<math.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<psio.h>
#include<libint.h>
#include<pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"

#include"dft_init.h"
#include"weighting.h"
#include"calc_den.h"
#include"functional.h"
#include"physconst.h"
#include"bragg.h"

void xc_fock(void){
  int i,j,k,l,m,n;
  int natoms, rpoints;
  int ang_points;
  int num_ao;

  int dum;
  int num;
  
  int nstri;
  double temp;
  double **tmpmat1;
  double *Gtri, *Gtri_o;  /* Total and open-shell 
			     G matrices and lower 
			     triagonal form in SO basis */  
  double r;
  double rind;
  double rpoints_double;
  double bragg;
  double qr;
  double drdq;
  double jacobian;
  struct coordinates geom;
  double xa,ya,za;
  double Becke_weight;
  double ang_quad;
  double four_pi_div_by_rps;
  double vfunc_val;
  double efunc_val;
  double vval;
  double eval = 0.0;
  double bas1,bas2;
    
  struct den_info_s den_info;

  natoms = Molecule.num_atoms;
  num_ao = BasisSet.num_ao;
  rpoints = DFT_options.rpoints;
  ang_points = DFT_options.grid_info[0].angpoints;
  rpoints_double = (double) rpoints;
  four_pi_div_by_rps = 4.0*_pi/rpoints_double;
  DFT_options.basis = (double *)malloc(sizeof(double)*num_ao);
  G = block_matrix(num_ao,num_ao);
  if(UserOptions.reftype == uhf)
      Go = block_matrix(num_ao,num_ao);
  
  /* Must set the Bragg-Slater radii */
  DFT_options.bragg = init_array(Molecule.num_atoms);
  
  for(i=0;i<Molecule.num_atoms;i++)
      DFT_options.bragg[i] =
	  Bragg_radii[(int) Molecule.centers[i].Z_nuc]*1.8897269;
  
  for(i=0;i<natoms;i++){
      xa = Molecule.centers[i].x;
      ya = Molecule.centers[i].y;
      za = Molecule.centers[i].z;
      
      bragg = DFT_options.bragg[i];
      
      for(j=1;j<rpoints;j++){
	  rind = (double) j;
	  
	  qr = j/ rpoints_double;
	  
	  /* -------------------------------
	     Straight from the Murray, Handy, Laming paper
	     for mr = 2 
	     ----------------------------------*/
	  
	  r = bragg*rind*rind/((rpoints_double - rind)
			       *(rpoints_double - rind));
	  
	  drdq = four_pi_div_by_rps*r*r*bragg*2.0*qr/((1-qr)*(1-qr)*(1-qr));
	  
	  
	  for(k=0;k<ang_points;k++){
	      
	      /* ----------------------------------
		 Calculate the cartesian points of the point
		 relative to the center of mass
		 -------------------------------------*/
      
	      geom.x = r*DFT_options.grid_info[0].leb_point[k].p_cart.x+xa;
	      geom.y = r*DFT_options.grid_info[0].leb_point[k].p_cart.y+ya;
	      geom.z = r*DFT_options.grid_info[0].leb_point[k].p_cart.z+za;
	      
	      
	      /*-----------------------------------
		Calculate the weighting funtion 
		----------------------------------*/
	      Becke_weight = weight_calc(i,geom,3);
	      
	      /*-----------------------------------
		Get the density information for this 
		point
		----------------------------------*/
	      
	      den_info = calc_density(geom);
	      
	      /*-------------------------------------
		Weight from Lebedev
		-----------------------------------*/
	      
	      ang_quad = DFT_options.grid_info[0].leb_point[k].ang_quad_weight;
	      
	      /*------------------------------------
		Calculate the potential functional
		and energy functional at this
		point
		-----------------------------------*/
	      
	      vfunc_val = DFT_options.exchange_V_function(den_info)
		  +DFT_options.correlation_V_function(den_info);
	      
	      efunc_val = DFT_options.exchange_function(den_info)
		  + DFT_options.correlation_function(den_info);
	      
	      vval = drdq*ang_quad*Becke_weight*vfunc_val;
	      eval += drdq*ang_quad*Becke_weight*efunc_val;
	      
	      /* ------------------------------------
		 Update the G matrix
		 -----------------------------------*/
	      
	      for(m=0;m<num_ao;m++){
		  bas1 = DFT_options.basis[m];
		  for(n=0;n<num_ao;n++){
		      bas2 = DFT_options.basis[n];
		      G[m][n] += vval*bas1*bas2;
		  }
	      }
	  }
      }
  }
  
  /*----------------------
    Transform to SO basis
    ----------------------*/
  if (Symmetry.nirreps > 1 || BasisSet.puream) {
      tmpmat1 = block_matrix(Symmetry.num_so,BasisSet.num_ao);
      mmult(Symmetry.usotao,0,G,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,G,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
      if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
	  mmult(Symmetry.usotao,0,Go,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
	  mmult(tmpmat1,0,Symmetry.usotao,1,Go,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
      }
      free_block(tmpmat1);
  }
  
  /*-------------------------
    Write G-matrices to disk
    -------------------------*/
  
  nstri = ioff[Symmetry.num_so];
  Gtri = init_array(nstri);
  sq_to_tri(G,Gtri,Symmetry.num_so);
  free_block(G);
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_write_entry(IOUnits.itapDSCF,"DFT XC-energy",(char *) &(eval), sizeof(double));
  switch (UserOptions.reftype) {
  case rohf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      psio_write_entry(IOUnits.itapDSCF, "Open-shell XC G-matrix"
		       , (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri_o);
  case rhf:
      psio_write_entry(IOUnits.itapDSCF, "Total XC G-matrix"
		       , (char *) Gtri, sizeof(double)*nstri);
      free(Gtri);
      break;
      
  case uhf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      /*--- Form alpha and beta Fock matrices first
	and then write them out ---*/      
      
      for(i=0;i<nstri;i++) {
	  temp = Gtri[i] + Gtri_o[i];
	  Gtri[i] = Gtri[i] - Gtri_o[i];
	  Gtri_o[i] = temp;
      }
      psio_write_entry(IOUnits.itapDSCF, "Alpha XC G-matrix"
		       , (char *) Gtri, sizeof(double)*nstri);
      psio_write_entry(IOUnits.itapDSCF, "Beta XC G-matrix"
		       , (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri);
      free(Gtri_o);
      break;
  }
  free(DFT_options.basis);
  psio_close(IOUnits.itapDSCF, 1);

  return;
}




