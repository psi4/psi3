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
#include"calc_den_fast.h"
#include"calc_den.h"
#include"functional.h"
#include"physconst.h"
#include"grid_init.h"
#include"free_grid_structs.h"
#include"dcr.h"
#include"init_close_shell_info.h"
#include"calc_close_basis.h"

 
void xc_fock(void){
  int i,j,k,l,m,n,q,s,t,u;
  int ua, atom, ua_deg;
  int rpoints;
  int ang_points;
  int num_ao;
  int point_count=0;
  int dum;
  int num;
  int moff,noff,mtmp,ntmp;
  int am2shell1,am2shell2;
  int shell_type1,shell_type2;
  int ioff1,ioff2;
  int chek = 1;
  int close_shells;
  int close_aos;
  int tmp;
  
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
  double ua_deg_d;
  double bragg;
  double qr;
  double drdq;
  double jacobian;
  double xa,ya,za;
  double Becke_weight;
  double ang_quad;
  double four_pi_div_by_rps;
  double exch_vfunc_val=0.0;
  double exch_efunc_val=0.0;
  double corr_vfunc_val=0.0;
  double corr_efunc_val=0.0;
  double den_val=0.0;
  double exch_vval=0.0;
  double corr_vval=0.0;
  double exch_eval=0.0;
  double corr_eval=0.0;
  double vval = 0.0;
  double eval = 0.0;
  double bas1 = 0.0;
  double bas2 = 0.0;
  double vvalbas = 0.0;
  
  struct coordinates geom;
  struct den_info_s den_info;
  
  struct atomic_grid_s *atm_grd;
  struct leb_chunk_s *chnk;
  leb_sphere_t *sphr;
  leb_point_t *pnt;
  
  
  fprintf(outfile,"\nPade_int = %10.10lf\n",Pade_int(0.5,-0.10498,3.72744,12.9352,0.06218414));
  num_ao = BasisSet.num_ao;
  DFT_options.basis = (double *)malloc(sizeof(double)*num_ao);
  
  G = init_matrix(num_ao,num_ao);
  if(UserOptions.reftype == uhf)
      Go = block_matrix(num_ao,num_ao);
  
  /* ----Initialize Close shell data structure ---- */
  
  DFT_options.close_shell_info = init_close_shell_info();

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
      ua_deg_d = (double) ua_deg;

      xa = atm_grd->atom_center.x;
      ya = atm_grd->atom_center.y;
      za = atm_grd->atom_center.z;
      
      bragg = atm_grd->Bragg_radii;
      
      for(j=0;j<atm_grd->chunk_num;j++){
	  chnk = &(atm_grd->leb_chunk[j]);
	  timer_on("close basis");
	  calc_close_basis(ua,j);
	  close_shells = DFT_options.close_shell_info.num_close_shells;
	  close_aos = DFT_options.close_shell_info.num_close_aos;
	  timer_off("close basis");
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
		      timer_on("DEN1");
		      den_info = calc_density_fast(geom,ua,j);
		      timer_off("DEN1");
		      /*timer_on("DEN2");
		      den_info = calc_density(geom);
		      timer_off("DEN2");*/
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
			  /*fprintf(outfile,"\nua_deg = %10.10lf",ua_deg_d);*/
			  den_val += 2.0*ua_deg_d*drdq*ang_quad
			      *Becke_weight*den_info.den;
			  
			  exch_vfunc_val = DFT_options.
			      exchange_V_function(den_info);
			      
			  corr_vfunc_val = DFT_options.
			      correlation_V_function(den_info);
		      
			  exch_efunc_val = DFT_options.
			      exchange_function(den_info);
			  corr_efunc_val = DFT_options.
			      correlation_function(den_info);
		      
			  exch_vval = ua_deg_d*drdq*ang_quad
			      *Becke_weight*exch_vfunc_val;
			  corr_vval = ua_deg_d*drdq*ang_quad
			      *Becke_weight*corr_vfunc_val;
			  vval = exch_vval+corr_vval;
			  
			  exch_eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*exch_efunc_val;
			  corr_eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*corr_efunc_val;
			  eval += ua_deg_d*drdq*ang_quad
			      *Becke_weight*(exch_efunc_val
					     +corr_efunc_val);
			  
			 
			  
			  /*------------------------------------
			    Update the G matrix
			    -----------------------------------*/
			  timer_on("FOCK");
			  t=0;
			  
			  for(m=0;m<close_aos;m++){
			      bas1 = vval*DFT_options.basis[m];
			      moff = DFT_options.close_shell_info.
				  aos_close_to_chunk[m];
			      for(n=m;n<close_aos;n++){
				  bas2 = bas1*DFT_options.basis[n];
				  noff = DFT_options.close_shell_info.
				      aos_close_to_chunk[n];
				  if(noff > moff)
				      G[moff][noff] += bas2;
				  else
				      G[noff][moff] += bas2;
			      }
			  }
			  timer_off("FOCK");
		      }

		  }
	      }
	  }
      }
  }

  
  free_close_shell_info(DFT_options.close_shell_info);
  /*print_mat(G,num_ao,num_ao,outfile);*/
   for(m=0;m<num_ao;m++)
      for(n=m;n<num_ao;n++)
	  G[n][m]=G[m][n];
  
  timer_off("DFT");
  timer_done();
  
  
  /*----------------------
    Unsort the Fock matrix
    back to shell ordering
    ---------------------*/
  
  
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
  fprintf(outfile,"\nX-Energy = %10.10lf",exch_eval);
  fprintf(outfile,"\nC-Energy = %10.10lf",corr_eval);
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
  cleanup_grid_type(DFT_options.grid);
  psio_close(IOUnits.itapDSCF, 1);
  return;
}




