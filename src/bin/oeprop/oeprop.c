#include "includes.h"
#include "oeprop.h"
#include "oeprop.gbl"


void print_intro();
void print_tasks();
void print_params();
void print_pop_header();
void print_mp();
void print_esp();
void print_grid();
void print_misc();
FILE *grid_file;

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
 wrt_dipints = 0;
 dip_file = 59;
 nuc_esp = 1;
 grid = 0;
 nix = 10;
 niy = 10;
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
 
 title = file30_rd_title();
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

 compute_onecgt();

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

 free_matrix(scf_evec_so,nbfso);
 free_matrix(usotao,nbfso);
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

void print_intro()
{ 
  fprintf(outfile,"    **********************************************\n");
  fprintf(outfile,"    *                    OEPROP                  *\n");
  fprintf(outfile,"    *          A simple property program         *\n");
  fprintf(outfile,"    *              by a big TOOL fan             *\n");
  fprintf(outfile,"    **********************************************\n\n");
}



void print_tasks()
{ 
   fprintf(outfile,"\n  TASKS to be performed :\n");
   
   if (read_opdm) {
     fprintf(outfile,"    $One-particle density in %s basis in %s form will be read from file%d",
             opdm_basis,opdm_format,opdm_file);
     if (asymm_opdm)
       fprintf(outfile," and symmetrized.\n");
     else
       fprintf(outfile,".\n");

     if (wrtnos)
       fprintf(outfile,"    $Natural orbitals will be written to file30.\n");
   }
   else
     fprintf(outfile,"    $One-particle density will be computed from the eigenvector in file30.\n");

   if (spin_prop)
     fprintf(outfile,"    $Spin properties will be evaluated.\n");
   
   switch (mpmax) {
     
     case 1:
       fprintf(outfile,"    $Only electric dipole moment will be computed.\n");
       break;
     
     case 2:
       fprintf(outfile,"    $Electric dipole and quadrupole moments will be computed.\n");
       break;
     
     case 3:
       fprintf(outfile,"    $Electric dipole, quadrupole, and octopole moments will be computed.\n");
       break;
   }
   
   fprintf(outfile,"    $Reference point for the electric multipole moments calculation is ");
   switch (mp_ref) {
     
     case 1: fprintf(outfile,"\n      the center of mass.\n");
             break;
             
     case 2: fprintf(outfile,"\n      the origin of the coordinate system.\n");
             break;

     case 3: fprintf(outfile,"\n      the center of electronic charge computed from a Mulliken analysis.\n");
             break;
     
     case 4: fprintf(outfile,"\n      the center of nuclear charge.\n");
             break;
     
     case 5: fprintf(outfile,"\n      the center of net charge computed from a Mulliken analysis.\n");
             break;
     
     default: fprintf(outfile,"\n      at (%lf %lf %lf)\n",
                      mp_ref_xyz[0],mp_ref_xyz[1],mp_ref_xyz[2]);
   }
  if (corr)
    fprintf(outfile,"    $Correlation corrections to electric multipole moments will be computed.\n");
  if (nuc_esp)
    fprintf(outfile,"    $Electrostatic properties at the nuclei will be evaluated.\n");
  if (grid) {
    switch(grid) {
      case 1: fprintf(outfile,"    $Electrostatic potential ");
  	      break;
      case 2: if (!spin_prop)
	        fprintf(outfile,"    $Electron density ");
              else
	        fprintf(outfile,"    $Spin density ");
	      break;
      case 3: if (!spin_prop)
	        fprintf(outfile,"    $Electron density gradient ");
              else
	        fprintf(outfile,"    $Spin density gradient ");
	      break;
      case 4: if (!spin_prop)
	        fprintf(outfile,"    $Laplacian of electron density ");
              else
	        fprintf(outfile,"    $Laplacian of spin density ");
	      break;
    }
    fprintf(outfile,"will be evaluated over a rectangular %dx%d grid.\n",nix+1,niy+1);
  }
  if (wrt_dipints)
    fprintf(outfile,"    $Dipole moment integrals will be written to file%d\n",dip_file);

}



void print_params()
{
   int i;

   fprintf(outfile,"\n  Title : '%s'\n",title);
   fprintf(outfile,"\n  List of PARAMETERS :\n");
   fprintf(outfile,"    # of atoms                 =\t%d\n",natom);
   fprintf(outfile,"    # of molecular orbitals    =\t%d\n",nmo);
   fprintf(outfile,"    # of basis functions       =\t%d\n",nbfso);
   fprintf(outfile,"    # of atomic orbitals       =\t%d\n",nbfao);
   fprintf(outfile,"    # of irreps                =\t%d\n",nirreps);
   fprintf(outfile,"    Total charge               =\t%d\n",charge);
   fprintf(outfile,"    # of unique shells         =\t%d\n",nshell);
   fprintf(outfile,"    # of primitives            =\t%d\n",nprim);
   fprintf(outfile,"    Print level                =\t%d\n",print_lvl);
   if (grid) {
     fprintf(outfile,"\n  List of GRID PARAMETERS :\n");
     fprintf(outfile,"    GRID_ORIGIN                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_origin[0],grid_origin[1],grid_origin[2]);
     fprintf(outfile,"    GRID_XY0                   =\t( %8.5lf %8.5lf )\n",grid_xy0[0],grid_xy0[1]);
     fprintf(outfile,"    GRID_XY1                   =\t( %8.5lf %8.5lf )\n",grid_xy1[0],grid_xy1[1]);
     fprintf(outfile,"    GRID_UNIT_X                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_unit_x[0],grid_unit_x[1],grid_unit_x[2]);
     fprintf(outfile,"    GRID_UNIT_Y                =\t( %8.5lf %8.5lf %8.5lf )\n",grid_unit_y[0],grid_unit_y[1],grid_unit_y[2]);
     fprintf(outfile,"    GRID_ZMIN                  =\t  %8.5lf\n",grid_zmin);
     fprintf(outfile,"    GRID_ZMAX                  =\t  %8.5lf\n",grid_zmax);
   }
   fprintf(outfile,"\n");
}


void print_pop_header()
{
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"   ** Mulliken population analysis of one-particle density **\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
}


void print_mp()
{
  FILE *file_dipmom;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"                *** Electric multipole moments ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  if (charge != 0) {
    fprintf(outfile,"  CAUTION : The system has non-vanishing charge, therefore dipole\n");
    fprintf(outfile,"    and higher moments depend on the reference point. \n\n");
  }
  else
  if ((dtot > 1.0E-15) && (mpmax > 1)) {
    fprintf(outfile,"  CAUTION : The system has non-vanishing dipole moment, therefore\n");
    fprintf(outfile,"    quadrupole and higher moments depend on the reference point.\n\n");
  }
  else
  if (((qvals[0]*qvals[0] + qvals[1]*qvals[1] + qvals[2]*qvals[2]) > 1.0E-15) 
      && (mpmax > 2)) {
    fprintf(outfile,"  CAUTION : The system has non-vanishing quadrupole moment, therefore\n");
    fprintf(outfile,"    octopole and higher moments depend on the reference point.\n\n");
  }

  fprintf(outfile," -Coordinates of the reference point (a.u.) :\n");
  fprintf(outfile,"           x                     y                     z\n");
  fprintf(outfile,"  --------------------  --------------------  --------------------\n");
  fprintf(outfile,"  %20.10lf  %20.10lf  %20.10lf\n\n",
          mp_ref_xyz[0],mp_ref_xyz[1],mp_ref_xyz[2]);
  if (print_lvl >= PRINTDIPCOMPLEVEL) {
    fprintf(outfile," -Contributions to electric dipole moment (a.u.) :\n\n");
    fprintf(outfile,"   -Electronic part :\n\n");
    fprintf(outfile,"    mu(X) =  %11.8lf  mu(Y) =  %11.8lf  mu(Z) = %11.8lf\n\n",
            dx_e,dy_e,dz_e);
    fprintf(outfile,"   -Nuclear part :\n\n");
    fprintf(outfile,"    mu(X) =  %11.8lf  mu(Y) =  %11.8lf  mu(Z) = %11.8lf\n\n",
            dx_n,dy_n,dz_n);
  }
  fprintf(outfile," -Electric dipole moment (expectation values) :\n\n");
  fprintf(outfile,"    mu(X)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dx*_dipmom_au2debye,dx*_dipmom_au2si,dx);
  fprintf(outfile,"    mu(Y)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dy*_dipmom_au2debye,dy*_dipmom_au2si,dy);
  fprintf(outfile,"    mu(Z)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dz*_dipmom_au2debye,dz*_dipmom_au2si,dz);
  fprintf(outfile,"    |mu|   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
          dtot*_dipmom_au2debye,dtot*_dipmom_au2si,dtot);
  if (corr) {
    fprintf(outfile,"\n -Correlation correction to electric dipole moment :\n\n");
    fprintf(outfile,"    cc(X)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dxcc*_dipmom_au2debye,dxcc*_dipmom_au2si,dxcc);
    fprintf(outfile,"    cc(Y)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dycc*_dipmom_au2debye,dycc*_dipmom_au2si,dycc);
    fprintf(outfile,"    cc(Z)  =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dzcc*_dipmom_au2debye,dzcc*_dipmom_au2si,dzcc);
    fprintf(outfile,"\n -Corrected electric dipole moment :\n\n");
    fprintf(outfile,"   mu(X)   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            (dx+dxcc)*_dipmom_au2debye,(dx+dxcc)*_dipmom_au2si,(dx+dxcc));
    fprintf(outfile,"   mu(Y)   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            (dy+dycc)*_dipmom_au2debye,(dy+dycc)*_dipmom_au2si,(dy+dycc));
    fprintf(outfile,"   mu(Z)   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            (dz+dzcc)*_dipmom_au2debye,(dz+dzcc)*_dipmom_au2si,(dz+dzcc));
    dtot = sqrt((dx+dxcc)*(dx+dxcc) + (dy+dycc)*(dy+dycc) + (dz+dzcc)*(dz+dzcc));
    fprintf(outfile,"    |mu|   =  %8.5lf D  =  %15.8e C*m  =  %11.8lf a.u.\n",
            dtot*_dipmom_au2debye,dtot*_dipmom_au2si,dtot);
  }

  /* write the total dipole moment to an ASCII file */
  if (wrt_dipmom) { 
    ffile(&file_dipmom,"dipmom.dat",1);
    fprintf(file_dipmom,"%20.10lf%20.10lf%20.10lf\n",
      (dx+dzcc)*_dipmom_au2debye,(dy+dycc)*_dipmom_au2debye,
      (dz+dzcc)*_dipmom_au2debye);
    fclose(file_dipmom);
  }


  if (mpmax > 1) {
    fprintf(outfile,"\n -Components of electric quadrupole moment (expectation values) (a.u.) :\n\n");
    fprintf(outfile,"     Q(XX) =  %12.8lf   Q(YY) =  %12.8lf   Q(ZZ) =  %12.8lf\n",
            qxx,qyy,qzz);
    fprintf(outfile,"     Q(XY) =  %12.8lf   Q(XZ) =  %12.8lf   Q(YZ) =  %12.8lf\n",
            qxy,qxz,qyz);
    if (corr) {
      fprintf(outfile,"\n -Correlation correction to electric quadrupole moment (a.u.) :\n\n");
      fprintf(outfile,"    cc(XX) =  %12.8lf  cc(YY) =  %12.8lf  cc(ZZ) =  %12.8lf\n",
              qxxcc,qyycc,qzzcc);
      fprintf(outfile,"    cc(XY) =  %12.8lf  cc(XZ) =  %12.8lf  cc(YZ) =  %12.8lf\n",
              qxycc,qxzcc,qyzcc);
      fprintf(outfile,"\n -Principal values (a.u.) and axis of corrected electric quadrupole moment :\n\n");
    }
    else
      fprintf(outfile,"\n -Principal values (a.u.) and axis of electric quadrupole moment :\n\n");
    fprintf(outfile,"    Q1     =  %12.8lf      V1 = (%11.8lf %11.8lf %11.8lf)\n",
            qvals[0],qvecs[0][0],qvecs[1][0],qvecs[2][0]);
    fprintf(outfile,"    Q2     =  %12.8lf      V2 = (%11.8lf %11.8lf %11.8lf)\n",
            qvals[1],qvecs[0][1],qvecs[1][1],qvecs[2][1]);
    fprintf(outfile,"    Q3     =  %12.8lf      V3 = (%11.8lf %11.8lf %11.8lf)\n",
            qvals[2],qvecs[0][2],qvecs[1][2],qvecs[2][2]);
  }
  if (mpmax > 2) {
    fprintf(outfile,"\n -Components of electric octopole moment (expectation values) (a.u.) :\n\n");
    fprintf(outfile,"    O(XXX) =  %12.8lf  O(XXY) =  %12.8lf  O(XXZ) =  %12.8lf\n",
            oxxx,oxxy,oxxz);
    fprintf(outfile,"    O(YYY) =  %12.8lf  O(XYY) =  %12.8lf  O(YYZ) =  %12.8lf\n",
            oyyy,oxyy,oyyz);
    fprintf(outfile,"    O(ZZZ) =  %12.8lf  O(XZZ) =  %12.8lf  O(YZZ) =  %12.8lf\n",
            ozzz,oxzz,oyzz);
    fprintf(outfile,"                            O(XYZ) =  %12.8lf \n",oxyz);
    if (corr) {
      fprintf(outfile,"\n -Correlation correction to electric octopole moment (a.u.) :\n\n");
      fprintf(outfile,"   cc(XXX) =  %12.8lf cc(XXY) =  %12.8lf cc(XXZ) =  %12.8lf\n",
              oxxxcc,oxxycc,oxxzcc);
      fprintf(outfile,"   cc(YYY) =  %12.8lf cc(XYY) =  %12.8lf cc(YYZ) =  %12.8lf\n",
              oyyycc,oxyycc,oyyzcc);
      fprintf(outfile,"   cc(ZZZ) =  %12.8lf cc(XZZ) =  %12.8lf cc(YZZ) =  %12.8lf\n",
              ozzzcc,oxzzcc,oyzzcc);
      fprintf(outfile,"                           cc(XYZ) =  %12.8lf \n",oxyzcc);
      fprintf(outfile,"\n -Corrected electric octopole moment (a.u.) :\n\n");
      fprintf(outfile,"    O(XXX) =  %12.8lf  O(XXY) =  %12.8lf  O(XXZ) =  %12.8lf\n",
              oxxx+oxxxcc,oxxy+oxxycc,oxxz+oxxzcc);
      fprintf(outfile,"    O(YYY) =  %12.8lf  O(XYY) =  %12.8lf  O(YYZ) =  %12.8lf\n",
              oyyy+oyyycc,oxyy+oxyycc,oyyz+oyyzcc);
      fprintf(outfile,"    O(ZZZ) =  %12.8lf  O(XZZ) =  %12.8lf  O(YZZ) =  %12.8lf\n",
              ozzz+ozzzcc,oxzz+oxzzcc,oyzz+oyzzcc);
      fprintf(outfile,"                            O(XYZ) =  %12.8lf \n",oxyz+oxyzcc);
    }
  }
  fprintf(outfile,"\n\n");
}


void print_esp()
{
  int i;
  
  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"      *** Electrostatic  properties at atomic centers ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of atomic centers (a.u.):\n");
  fprintf(outfile,"    #   Charge           x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ---  ------  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%5d%7d    %20.10lf  %20.10lf  %20.10lf\n",
                    i+1,(int)zvals[i],geom[i][0]+mp_ref_xyz[0],
                        geom[i][1]+mp_ref_xyz[1],
                        geom[i][2]+mp_ref_xyz[2]);
  fprintf(outfile,"\n\n");
  fprintf(outfile," -Electrostatic potential and electric field (a.u.) :\n\n");
  fprintf(outfile,"    Center         phi            Ex             Ey           Ez\n");
  fprintf(outfile,"    ------    ------------   ------------  ------------  ------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d      %12.8lf   %12.8lf  %12.8lf  %12.8lf\n",
            i+1,phi[i],ex[i],ey[i],ez[i]);
  fprintf(outfile,"\n\n");
  fprintf(outfile," -Electric field gradient (regular form) (a.u.):\n\n");
  fprintf(outfile,"    Center           XX                    YY                    ZZ\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
            i+1,dexx[i],deyy[i],dezz[i]);
  fprintf(outfile,"\n    Center           XY                    XZ                    YZ\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
            i+1,dexy[i],dexz[i],deyz[i]);
  fprintf(outfile,"\n\n");
  fprintf(outfile," -Electric field gradient (traceless tensor form) (a.u.):\n\n");
  fprintf(outfile,"    Center        XX - RR/3             YY - RR/3             ZZ - RR/3\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",i+1,
            (2*dexx[i]-deyy[i]-dezz[i])/3,
            (2*deyy[i]-dexx[i]-dezz[i])/3,
            (2*dezz[i]-dexx[i]-deyy[i])/3);
  fprintf(outfile,"\n    Center           XY                    XZ                    YZ\n");
  fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
            i+1,dexy[i],dexz[i],deyz[i]);
  fprintf(outfile,"\n\n");

  if (spin_prop) {
    fprintf(outfile," -Dipole-dipole term in hyperfine coupling constant (tensor form) (a.u.):\n\n");
    fprintf(outfile,"    Center        XX - RR/3             YY - RR/3             ZZ - RR/3\n");
    fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
              i+1,(2*ahfsxx[i]-ahfsyy[i]-ahfszz[i])/3,(2*ahfsyy[i]-ahfsxx[i]-ahfszz[i])/3,(2*ahfszz[i]-ahfsxx[i]-ahfsyy[i])/3);
    fprintf(outfile,"\n    Center           XY                    XZ                    YZ\n");
    fprintf(outfile,"    ------  --------------------  --------------------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf  %20.8lf  %20.8lf\n",
	      i+1,ahfsxy[i],ahfsxz[i],ahfsyz[i]);
    fprintf(outfile,"\n\n");
  }

  if (spin_prop) {
    fprintf(outfile," -Electron and spin densities (a.u.):\n\n");
    fprintf(outfile,"    Center    Electron density         Spin density\n");
    fprintf(outfile,"    ------  --------------------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf  %20.8lf\n",i+1,edens[i],sdens[i]);
  }
  else {
    fprintf(outfile," -Electron density (a.u.):\n\n");
    fprintf(outfile,"    Center           rho\n");
    fprintf(outfile,"    ------  --------------------\n");
    for(i=0;i<natom;i++)
      fprintf(outfile,"%8d    %20.8lf\n",i+1,edens[i]);
  }
  fprintf(outfile,"\n\n");
}


void print_grid()
{
  int i,j,k;
  double step_x, step_y, x, y;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"    *** Evaluating properties over a rectangular 2D grid ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");
  fprintf(outfile," -Coordinates of the lower left, lower right, and upper left corners of\n");
  fprintf(outfile,"  the grid rectangle (a.u.):\n");
  fprintf(outfile,"    **            x");
  fprintf(outfile,"                     y                     z\n");
  fprintf(outfile,"   ----  --------------------");
  fprintf(outfile,"  --------------------  --------------------\n");
  fprintf(outfile,"    LL   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0],grid_origin[1],grid_origin[2]);
  fprintf(outfile,"    LR   %20.10lf  %20.10lf  %20.10lf\n",
	  grid_origin[0]+grid_unit_x[0]*(grid_xy1[0]-grid_xy0[0]),
	  grid_origin[1]+grid_unit_x[1]*(grid_xy1[0]-grid_xy0[0]),
	  grid_origin[2]+grid_unit_x[2]*(grid_xy1[0]-grid_xy0[0]));
  fprintf(outfile,"    UL   %20.10lf  %20.10lf  %20.10lf\n\n\n",
	  grid_origin[0]+grid_unit_y[0]*(grid_xy1[1]-grid_xy0[1]),
	  grid_origin[1]+grid_unit_y[1]*(grid_xy1[1]-grid_xy0[1]),
	  grid_origin[2]+grid_unit_y[2]*(grid_xy1[1]-grid_xy0[1]));
  
  switch (grid) {
    case 1:
      grid_file = fopen("esp.dat","w");
      break;

    case 2:
      if (!spin_prop)
        grid_file = fopen("edens.dat","w");
      else
	grid_file = fopen("sdens.dat","w");
      break;
    
    case 3:
      if (!spin_prop)
        grid_file = fopen("edgrad.dat","w");
      else
	grid_file = fopen("sdgrad.dat","w");
      break;

    case 4:
      if (!spin_prop)
        grid_file = fopen("edlapl.dat","w");
      else
	grid_file = fopen("sdlapl.dat","w");
      break;
  }

/*    fprintf(grid_file,"%8.4lf %d %8.4lf  %8.4lf %d %8.4lf  %8.4lf %d %8.4lf\n",
            grid_xyz0[0],nix,grid_xyz1[0],
            grid_xyz0[1],niy,grid_xyz1[1],
            grid_xyz0[2],niz,grid_xyz1[2]); */
  switch (grid) {
    case 1:
    case 2:
    case 4:
      fprintf(grid_file,"$DATA = CONTOUR\n");
      fprintf(grid_file,"%% xmin = %lf xmax = %lf nx = %d\n",grid_xy0[0],grid_xy1[0],nix+1);
      fprintf(grid_file,"%% ymin = %lf ymax = %lf ny = %d\n",grid_xy0[1],grid_xy1[1],niy+1);
      fprintf(grid_file,"%% zmin = %lf zmax = %lf\n",grid_zmin,grid_zmax);
      fprintf(grid_file,"%% contfill = T meshplot = T\n");
      for(i=0;i<=niy;i++) {
        for(j=0;j<=nix;j++)
	  if (grid_pts[j][i] < grid_zmin)
            fprintf(grid_file," %lf ",grid_zmin);
	  else
	    if (grid_pts[j][i] > grid_zmax)
	      fprintf(grid_file," %lf ",grid_zmax);
	    else
	      fprintf(grid_file," %lf ",grid_pts[j][i]);
        fprintf(grid_file,"\n");
      }
      break; 
      
    case 3:
      fprintf(grid_file,"$DATA = VECTOR\n");
      fprintf(grid_file,"%% xmin = %lf xmax = %lf\n",grid_xy0[0],grid_xy1[0]);
      fprintf(grid_file,"%% ymin = %lf ymax = %lf\n",grid_xy0[1],grid_xy1[1]);
      fprintf(grid_file,"%% zmin = %lf zmax = %lf\n",0.0,0.0);
      fprintf(grid_file,"%% linecolor = 3\n");
      fprintf(grid_file,"%% xlog = off vscale = 0.2\n\n");
      step_x = (grid_xy1[0]-grid_xy0[0])/nix;
      step_y = (grid_xy1[1]-grid_xy0[1])/niy;
      for(i=0;i<=nix;i++) {
	x = grid_xy0[0] + step_x*i;
	for(j=0;j<=niy;j++) {
        y = grid_xy0[1] + step_y*j;
        if (fabs(grid_pts[i][j]) <= MAXDENSGRAD)
          fprintf(grid_file,"%9.5lf  %9.5lf  %9.5lf  %12.8lf  %12.8lf  %12.8lf\n",x,y,0.0,
                  grid_vecX[i][j],grid_vecY[i][j],0.0);
        }
      }
      break;
  }
      
  fprintf(grid_file,"$END\n");
  fclose(grid_file);
}


void print_misc()
{
  int i,j,k;

  fprintf(outfile," --------------------------------------------------------------\n");
  fprintf(outfile,"                *** Miscellaneous properties ***\n");
  fprintf(outfile," --------------------------------------------------------------\n\n");

  fprintf(outfile," -Relativistic MVD one-electron corrections to the energy (a.u.):\n\n");
  fprintf(outfile,"    Mass-velocity (p^4) term     :   %12.8lf\n",massveloc);
  fprintf(outfile,"    One-electron Darwin term     :   %12.8lf\n",darw);
  fprintf(outfile,"    Total one-electron MVD terms :   %12.8lf\n",massveloc+darw);
  fprintf(outfile,"\n\n");

  if (mpmax > 1) {
    fprintf(outfile,"  NOTE : Spatial extents are computed with respect to the same reference point\n");
    fprintf(outfile,"         as multipole moments.\n\n");
    fprintf(outfile," -Electronic spatial extents (a.u.) :\n\n");
    fprintf(outfile,"     <X^2> =  %11.4lf    <Y^2> =  %11.4lf    <Z^2> =  %11.4lf\n",
            exp_x2,exp_y2,exp_z2);
    fprintf(outfile,"                             <R^2> =  %11.4lf\n",
            exp_x2+exp_y2+exp_z2);
    fprintf(outfile,"\n -Orbital spatial extents ");
    if (read_opdm) 
      fprintf(outfile,"of NOs contructed from onepdm in file%d ",opdm_file);
    else
      fprintf(outfile,"of MOs in file30 ");
    fprintf(outfile,"(a.u.) :\n\n");
    fprintf(outfile,"    MO #   Symm     <X^2>        <Y^2>        <Z^2>        <R^2>\n");
    fprintf(outfile,"   ------  ----  -----------  -----------  -----------  -----------\n");
    k = 0;
    for(i=0;i<nirreps;i++)
      for(j=0;j<orbspi[i];j++)
        fprintf(outfile,"   %4d     %3s   %9.4lf    %9.4lf    %9.4lf    %9.4lf\n",
                k+1,irr_labs[i],MOXX[k],MOYY[k],MOZZ[k],MOXX[k]+MOYY[k]+MOZZ[k++]);
    fprintf(outfile,"\n");
  }
}


char *gprgid()
{
 char *prgid = "OEPROP";
   
 return(prgid);
}
