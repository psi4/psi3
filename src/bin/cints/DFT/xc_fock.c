#include<math.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<psio.h>
#include<libint.h>
#include<pthread.h>
#include<qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"

#include"dft_init.h"
#include"weighting.h"
#include"calc_den.h"
#include"functional.h"
#include"physconst.h"
#include"grid_init.h"
#include"dcr.h"

void xc_fock(void){
  int i,j,k,l,m,n;
  int ua, atom, ua_deg;
  int rpoints;
  int ang_points;
  int num_ao;
  int point_count=0;
  int dum;
  int num;
  
  int nstri;
  double temp;
  double *temp_arr;
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
  double xa,ya,za;
  double Becke_weight;
  double ang_quad;
  double four_pi_div_by_rps;
  double exch_vfunc_val;
  double exch_efunc_val;
  double corr_vfunc_val;
  double corr_efunc_val;
  double den_val;
  double exch_vval;
  double corr_vval;
  double exch_eval;
  double corr_eval;
  double vval = 0.0;
  double eval = 0.0;
  double bas1,bas2;
   
  struct coordinates geom;
  struct den_info_s den_info;
  
  struct atomic_grid_s *atm_grd;
  struct leb_chunk_s *chnk;
  leb_sphere_t *sphr;
  leb_point_t *pnt;

  num_ao = BasisSet.num_ao;
  DFT_options.basis = (double *)malloc(sizeof(double)*num_ao);
  
  G = block_matrix(num_ao,num_ao);
  if(UserOptions.reftype == uhf)
      Go = block_matrix(num_ao,num_ao);
  timer_init();
  timer_on("DFT");

  grid_init();
  
  /*-------------------------------------------------------
    Loop over symmetry-unique atoms only since integration
    domains around symm.-unique atoms are equivalent
    We are NOT employing the symmetry of angular grids
    about atoms in special positions (on symm. elements)
    like Handy & co. do
   -------------------------------------------------------*/
  for(ua=0;ua<Symmetry.num_unique_atoms;ua++){
      atm_grd = &(DFT_options.grid.atomic_grid[ua]);
      
      atom = atm_grd->atom_num;
      
/*--- Cheap trick to get degeneracies of each unique atom ---*/
      
      ua_deg = atm_grd->atom_degen;
     

      xa = atm_grd->atom_center.x;
      ya = atm_grd->atom_center.y;
      za = atm_grd->atom_center.z;
      
      bragg = atm_grd->Bragg_radii/*1.8897269*/;
      fprintf(outfile,"Bragg = %10.10lf",bragg);
      for(j=0;j<atm_grd->chunk_num;j++){
	  chnk = &(atm_grd->leb_chunk[j]);
	  
	  for(k=0;k<chnk->size;k++){
	      sphr = &(chnk->spheres[k]);
	  
	      r = sphr->r*bragg;
	      drdq = sphr->drdq*bragg*bragg*bragg;

	      for(l=0;l<sphr->n_ang_points;l++){
		  pnt = &(sphr->points[l]);
		  /* ----------------------------------
		     Calculate the cartesian points of the point
		     relative to the center of mass
		     -------------------------------------*/
      
		  geom.x = bragg*pnt->p_cart.x+xa;
		  geom.y = bragg*pnt->p_cart.y+ya;
		  geom.z = bragg*pnt->p_cart.z+za;
	      
		  /*-----------------------------------
		    Calculate the weighting funtion 
		    ----------------------------------*/
		  Becke_weight = weight_calc(atom,geom,3);
		  if(Becke_weight> WEIGHT_CUTOFF){
		      /*-----------------------------------
			Get the density information for this 
			point
			----------------------------------*/
		      den_info = calc_density(geom);
		      if(den_info.den > DEN_CUTOFF){
			  /*-------------------------------------
			    Weight from Lebedev
			    -----------------------------------*/
		      
			  ang_quad = pnt->ang_weight;
		      
			  /*------------------------------------
			    Calculate the potential functional
			    and energy functional at this
			    point
			    -----------------------------------*/
		      
			  den_val += 2.0*ua_deg*drdq*ang_quad
			      *Becke_weight*den_info.den;
			  
			  exch_vfunc_val = DFT_options.
			      exchange_V_function(den_info);
			      
			  corr_vfunc_val = DFT_options.
			      correlation_V_function(den_info);
		      
			  exch_efunc_val = DFT_options.
			      exchange_function(den_info);
			  corr_efunc_val = DFT_options.
			      correlation_function(den_info);
		      
			  exch_vval = ua_deg*drdq*ang_quad
			      *Becke_weight*exch_vfunc_val;
			  corr_vval = ua_deg*drdq*ang_quad
			      *Becke_weight*corr_vfunc_val;
			  vval = exch_vval+corr_vval;
			  
			  exch_eval += ua_deg*drdq*ang_quad
			      *Becke_weight*exch_efunc_val;
			  corr_eval += ua_deg*drdq*ang_quad
			      *Becke_weight*corr_efunc_val;
			  eval += ua_deg*drdq*ang_quad
			      *Becke_weight*(exch_efunc_val
					     +corr_efunc_val);
			  
			 
			  
			  /*------------------------------------
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
	  }
      }
  }
  timer_off("DFT");
  timer_done();
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
  fprintf(outfile,"\nDFT_energy = %10.10lf",eval);
  fprintf(outfile,"\ntrace of density = %10.10lf",den_val);
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_write_entry(IOUnits.itapDSCF,"DFT X-energy",
		   (char *) &(exch_eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT C-energy",
		   (char *) &(corr_eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT XC-energy",
		   (char *) &(eval), sizeof(double));
  psio_write_entry(IOUnits.itapDSCF,"DFT Den",
		   (char *) &(den_val), sizeof(double));
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
  /*-- Cleanup the DFT stuff and close files --*/
  /*cleanup_dft_options(DFT_options);*/
  psio_close(IOUnits.itapDSCF, 1);
  return;
}




