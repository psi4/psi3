#include "includes.h"
#include "prototypes.h"
#include "globals.h"


int main(int argc, char* argv) {

 int i,j,k,l,count;
 char buffer[80];	/* buffer string */

	/* Setting defaults */
 read_opdm = 0;
 opdm_file = 76;
 asymm_opdm = 0;
 wrtnos = 0;
 spin_prop = 0;
 print_lvl = 1;
 wrt_dipmom = 1;
 corr = 0;
 mpmax = 1;
 mp_ref = 0;
 nuc_esp = 1;
 grid = 0;
 grid3d = 0;
 nix = 10;
 niy = 10;
 niz = 10;
 grid_zmin = 0.0;
 grid_zmax = 3.0;
 edgrad_logscale = 5;
 zvec_file = 86;
 delete_zvec = 1;
 

	/* Initialization and printing intro */
 
 ffile(&infile,"input.dat",2);
 ffile(&outfile,"output.dat",1);
 ip_set_uppercase(1);
 ip_initialize(infile,outfile);
 ip_cwk_add(":DEFAULT");
 ip_cwk_add(":OEPROP");
 tstart(outfile);
 file30_init();
 print_intro();


/*************************** Main Code *******************************/

	/* Reading in basic information from file30 */
 
 title = file30_rd_label();
 natom = file30_rd_natom();
 natom3 = natom * 3;	/* wrong if the space is other than 3-dimensional */
 nmo = file30_rd_nmo();
 nbfso = file30_rd_nso();
 nbfao = file30_rd_nao();
 natri = nbfao * (nbfao+1)/2;
 nstri = nbfso * (nbfso+1)/2;
 nshell = file30_rd_nshell();
 nprim = file30_rd_nprim();
 iopen = file30_rd_iopen();
 nirreps = file30_rd_nirreps(); 
 nsym = file30_rd_nsymhf();
 orbspi = file30_rd_orbspi();
 clsdpi = file30_rd_clsdpi();    
 openpi = file30_rd_openpi();
 irr_labs = file30_rd_irr_labs();
 geom = file30_rd_geom();
 zvals = file30_rd_zvals();
 scf_evec_so = file30_rd_scf();
 usotao = file30_rd_usotao_new();
 scf_evec_ao = init_matrix(nbfao,nmo);
 mmult(usotao,1,scf_evec_so,0,scf_evec_ao,0,nbfao,nbfso,nmo,0);
    
	/* Parsing */

 parsing();


        /* Computing unit vectors for the 2D grid if neccessary */

 if (grid)
   grid_unitvec();


	/* Computing total charge of the system */

 charge = 0;
 for(i=0;i<nirreps;i++)
   charge -= 2*clsdpi[i] + openpi[i];
 for(i=0;i<natom;i++)
   charge += zvals[i];



	/* Setting up an offset array */

 ioff = init_int_array(nbfao);
 ioff[0] = 0;
 for(i=1;i<nbfao;i++) {
   ioff[i] = ioff[i-1] + i;
 }

	/* Computing double factorials df[i] = (i-1)!! */
 df[0] = 1.0;
 df[1] = 1.0;
 df[2] = 1.0;
 for(i=3; i<MAXFACT*2; i++) {
   df[i] = (i-1)*df[i-2];
 }
                 
                      

	/* Printing tasks and parameters */

 if (print_lvl >= PRINTTASKPARAMLEVEL) {
   print_tasks();
   print_params();
 }

	/* Reading in basis set inforamtion */

 read_basset_info();
 init_xyz(); 


	/* Obtain a density matrix */

 if (read_opdm)
   read_density();
 else
   compute_density();


	/* Obtain natural orbitals */

 if (read_opdm) 
   get_nmo(); 

 file30_close();

	/* Reading in Z-vector if neccessary */

 if (corr)
   read_zvec();


	/* Computing overlap matrix */

 compute_overlap();


	/* Mulliken population analysis */

 print_pop_header();
 populate();
 
	/* Compute coordinates of the MP reference point if needed */

 if (mp_ref != -1)
   compute_mp_ref_xyz();


	/* Moving molecule to the reference center!
	   Attention, ever since coordinates of atoms and of the grid box 
	   are stored in this new coordinate system. All coordinates are 
	   transformed back at the moment of printing out. */

 for(i=0;i<natom;i++) {
   geom[i][0] -= mp_ref_xyz[0];
   geom[i][1] -= mp_ref_xyz[1];
   geom[i][2] -= mp_ref_xyz[2];
 }
 if (grid) {
   grid_origin[0] -= mp_ref_xyz[0];
   grid_origin[1] -= mp_ref_xyz[1];
   grid_origin[2] -= mp_ref_xyz[2];
 }
 
	/* Computing one-electron integrals in 
	   terms of Cartesian Gaussians, 
	   electric first (dipole), second and third 
	   moments W.R.T. origin, electrostatic potential,
	   electric field and field gradients,
	   electron and spin density at atomic centers,
	   and various properties over a grid, if neccessary. */

 compute_oeprops();
 if (grid)
     compute_grid();

 print_mp();
 if (nuc_esp)
   print_esp();
 if (grid) {
   grid_origin[0] += mp_ref_xyz[0];
   grid_origin[1] += mp_ref_xyz[1];
   grid_origin[2] += mp_ref_xyz[2];
   print_grid();
 }
 print_misc();

	/* Cleaning up */

/* TDC --- converted libfile30 to use block_matrix() */
 free_block(scf_evec_so);
 free_block(usotao);
/* free_matrix(scf_evec_so,nbfso);
 free_matrix(usotao,nbfso); */
 free_matrix(scf_evec_ao,nbfao);
 free(ioff);
 free(Ptot);
 free(phi);
 free(ex); free(ey); free(ez);
 free(dexx); free(deyy); free(dezz);
 free(dexy); free(dexz); free(deyz);
 if (spin_prop) {
   free(Pspin);
   free(ahfsxx);
   free(ahfsyy);
   free(ahfszz);
   free(ahfsxy);
   free(ahfsxz);
   free(ahfsyz);
 }
 free(S);
 tstop(outfile);
 ip_done();
 exit(0);
 
}


char *gprgid()
{
 char *prgid = "OEPROP";
   
 return(prgid);
}
