#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <file30_params.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"


void print_intro();
void print_options();
void print_geometry(double);
void print_full_geometry(double);
void print_unique_geometry(double);
void print_symm_trans();
void print_basis_info();
void cleanup();

int main(int argc, char *argv[])
{
   /*variables and arrays*/
   int i,j,k,l;
   int errcod;
   double temp = 0.0;			/*temporary anything*/
   int natomtri = 0;                    /*number of elements in the lower triangle of a (num_atoms)x(num_atoms) matrix */
   double Z = 0.0;			/*Atomic Number*/
   double *Distance;			/*Lower triangle of the internuclear distance matrix*/
   double repulsion;			/*nuclear Repulsion Energy*/
   double ***E;				/*Unit Vectors between atoms*/
   double ***Bond_Angle;		/*Matrix holds bond angles*/
   FILE *pbasis = NULL;
   FILE *user_basis = NULL;
   FILE *local_basis = NULL;
   char *user_basis_filename, *user_basis_file;


     /*-------------------------------------
       Initialize files and parsing library
      -------------------------------------*/
     parsing_cmdline(argc,argv);
     start_io();
     pbasis = fopen(SHARE, "r");
     errcod = ip_string("BASIS_FILE",&user_basis_file,0);
     if (errcod == IPE_OK && strlen(user_basis_file) != 0) {
       /* Checking if only path to the file has been specified */
       if (!strncmp(user_basis_file+strlen(user_basis_file)-1,"/",1)) { /* Append default name - basis.dat */
	 user_basis_filename = malloc(sizeof(char)*(strlen(user_basis_file)+10));
	 strcpy(user_basis_filename,user_basis_file);
	 user_basis_filename = strcat(user_basis_filename,"basis.dat");
       }
       else
	 user_basis_filename = user_basis_file;
       user_basis = fopen(user_basis_filename, "r");
       free(user_basis_filename);
     }
     local_basis = fopen("./basis.dat", "r");
     if (local_basis != NULL) {
       ip_append(local_basis, outfile);
       fclose(local_basis);
     }
     if (user_basis != NULL) {
       ip_append(user_basis, outfile);
       fclose(user_basis);
     }
     if (pbasis != NULL) {
       ip_append(pbasis, outfile);
       fclose(pbasis);
     }
     ip_cwk_add(":BASIS");

     init_globals();
     parsing();
     print_intro();
     print_options();
     
     /*-----------------
       Read in geometry
      -----------------*/
     if (chkpt_geom == 0) { /* read geometry from input.dat */
	 if (cartOn)
	     read_cart();
	 else
	     read_zmat();
     }
     else { /* else read the next molecular geometry from checkpoint file */
	 read_chkpt_geom();
     }

     repulsion = 0.0;
     if (num_atoms > 1) {
       natomtri = num_atoms*(num_atoms+1)/2;
       Distance = init_array(natomtri);
       calc_distance(geometry,Distance,num_atoms);
       Nuc_repulsion(Distance, &repulsion);
     }

     fprintf(outfile,"\n  -Geometry before Center-of-Mass shift (a.u.):\n");
     print_geometry(1.0);

     reorient();
     fprintf(outfile,"\n  -Geometry after Center-of-Mass shift and reorientation (a.u.):\n");
     print_geometry(1.0);
     
     /*--------------
       Get symmetric
      --------------*/
     find_symmetry();
     count_uniques();
     print_symm_trans();

     /*---------------------------------------
       Canonical frame is the reference frame
       unless read geometry from chkpt file
      ---------------------------------------*/
     if (chkpt_geom == 0)
	 canon_eq_ref_frame();
     
     /*---------------------
       Parse basis set data
      ---------------------*/
     
     atom_basis = (char **) malloc(sizeof(char *)*num_atoms);
     for(i=0;i<num_atoms;i++)
       atom_basis[i] = NULL;
     read_basis();
     
     /*----------------------------------------------
       Form symmetry information arrays of all kinds
      ----------------------------------------------*/

     build_transmat();
     if (puream)
       build_cart2pureang();
     build_so_classes();
     build_usotao();

     /*-------------
       Read old MOs
      -------------*/
     if (chkpt_mos) {
	 init_oldcalc();
     }

     /*-------------------------------------------------
       Write the information out to the checkpoint file
      -------------------------------------------------*/
     write_to_file30(repulsion);

     /*-------------------------------------------------
       Project old MOs onto new basis and write 'em out
      -------------------------------------------------*/
     if (chkpt_mos) {
       oldcalc_projection();
       write_scf_to_file30();
       cleanup_oldcalc();
     }
     
     /*----------------------
       Print geometries etc.
      ----------------------*/

     print_basis_info();
     fprintf(outfile,"\n  -Unique atoms in the canonical coordinate system (a.u.):\n");
     print_unique_geometry(1.0);
     fprintf(outfile,"\n  -Geometry in the canonical coordinate system (a.u.):\n");
     print_geometry(1.0);
     fprintf(outfile,"\n  -Geometry in the canonical coordinate system (Angstrom):\n");
     print_geometry(_bohr2angstroms);
     rotate_geometry(geometry,Rref);
     fprintf(outfile,"\n  -Geometry in the reference coordinate system (a.u.):\n");
     print_geometry(1.0);
     if(!cartOn) {
	fprintf(outfile,"\n  -Full geometry in the canonical coordinate system (a.u.):\n");
	print_full_geometry(1.0);
     }

     if (num_atoms > 1) {
       /* Print the Nuclear Repulsion Energy */
       fprintf(outfile,"\n    Nuclear Repulsion Energy (a.u.) = %20.12lf\n\n",repulsion);

       /* Print the interatomic Distance Matrix in angstroms */
       for(i=0;i<natomtri;i++)
         Distance[i] *= _bohr2angstroms;
       fprintf(outfile,"  -The Interatomic Distances in angstroms:\n");
       print_array(Distance,num_atoms,outfile);
       free(Distance);
       fprintf(outfile,"\n\n");
     }


     cleanup();
     stop_io();
     exit(0);

}



void print_intro()
{
  if (print_lvl) {
      fprintf(outfile,"                                --------------\n");
      fprintf(outfile,"                                  WELCOME TO\n");
      fprintf(outfile,"                                    PSI  3\n");
      fprintf(outfile,"                                --------------\n\n");
      fflush(outfile);
  }

  return;
}


void print_options()
{
  if (print_lvl) {
      if (strlen(label))
	  fprintf(outfile,"  LABEL       =\t%s\n", label);
      fprintf(outfile,"  SHOWNORM    =\t%d\n",shownorm);
      fprintf(outfile,"  PUREAM      =\t%d\n",puream);
      if (subgroup != NULL) {
	  fprintf(outfile,"  SUBGROUP    =\t%s\n",subgroup);
	  if (unique_axis != NULL)
	      fprintf(outfile,"  UNIQUE_AXIS =\t%s\n",unique_axis);
      }
      fprintf(outfile,"  PRINT_LVL   =\t%d\n",print_lvl);
      if (chkpt_mos) {
	  fprintf(outfile,"  Will read old MO vector from the checkpoint file\n");
	  fprintf(outfile,"  and project onto the new basis.\n");
      }
      fflush(outfile);
  }

  return;
}


void print_geometry(double conv_factor)
{
  int i,j;

  fprintf(outfile,"       Center              X                  Y                   Z\n");
  fprintf(outfile,"    ------------   -----------------  -----------------  -----------------\n");

  for(i=0;i<num_atoms;i++){
    fprintf(outfile,"    %12s ",element[i]); fflush(outfile);
    for(j=0;j<3;j++)
      fprintf(outfile,"  %17.12lf",geometry[i][j]*conv_factor);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  fflush(outfile);

  return;
}


void print_full_geometry(double conv_factor)
{
  int i,j;

  fprintf(outfile,"     Center      Z           X                  Y                   Z\n");
  fprintf(outfile,"  ------------  ---  -----------------  -----------------  -----------------\n");

  for(i=0;i<num_entries;i++){
    fprintf(outfile,"  %12s %3d",full_element[i], (int) nuclear_charges[i]); fflush(outfile);
    for(j=0;j<3;j++)
      fprintf(outfile,"  %17.12lf",full_geom[i][j]*conv_factor);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"  --------------------------------------------------------------------------\n");
  fflush(outfile);

  return;
}


void print_unique_geometry(double conv_factor)
{
  int i,j;
  
  fprintf(outfile,"       Center              X                  Y                   Z\n");
  fprintf(outfile,"    ------------   -----------------  -----------------  -----------------\n");
  for(i=0;i<num_uniques;i++){
    fprintf(outfile,"    %12s ",element[u2a[i]]);
    for(j=0;j<3;j++)
      fprintf(outfile,"  %17.12lf",geometry[u2a[i]][j]*conv_factor);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  fflush(outfile);

  return;
}


void print_symm_trans()
{
  int i,j;

  fprintf(outfile,"\n  -SYMMETRY INFORMATION:\n");
  fprintf(outfile,"    Computational point group is %s\n",symmetry);
  fprintf(outfile,"    Number of irr. rep.      = %d\n",nirreps);
  fprintf(outfile,"    Number of atoms          = %d\n",num_atoms);
  fprintf(outfile,"    Number of unique atoms   = %d\n\n",num_uniques);
  if (print_lvl >= PRINTSYMINFO) {
    for(i=0;i<num_uniques;i++)
     fprintf(outfile,"    Unique atom %d generates %d atoms\n",i+1,unique_degen[i]);
    fprintf(outfile,"\n    Orbits of unique atoms:\n");
    for(i=0;i<num_uniques;i++) {
      fprintf(outfile,"     %d    ",u2a[i]+1);
      for(j=0;j<nirreps;j++)
	fprintf(outfile,"%d  ",atom_orbit[u2a[i]][j]+1);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    fprintf(outfile,"    Number of classes        = %d\n",num_classes);
    fprintf(outfile,"    Number of unique classes = %d\n\n",num_unique_classes);
    fprintf(outfile,"    Atom   Class   Symmetry position\n");
    fprintf(outfile,"    ----   -----   -----------------\n");
    for(i=0;i<num_atoms;i++)
      fprintf(outfile,"    %3d    %3d     %10d\n",i+1,atom_class[i]+1,atom_position[i]);
    fprintf(outfile,"\n");
    fprintf(outfile,"    Orbits of classes:\n");
    for(i=0;i<num_classes;i++) {
      fprintf(outfile,"     %d    ",i+1);
      for(j=0;j<nirreps;j++)
	fprintf(outfile,"%d  ",class_orbit[i][j]+1);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
  }
  fflush(outfile);

  return;
}


void print_basis_info()
{
  int i,j;

  fprintf(outfile,"  -BASIS SET INFORMATION:\n");
  fprintf(outfile,"    Total number of shells = %d\n",num_shells);
  fprintf(outfile,"    Number of primitives   = %d\n",num_prims);
  fprintf(outfile,"    Number of AO           = %d\n",num_ao);
  fprintf(outfile,"    Number of SO           = %d\n\n",num_so);
  fprintf(outfile,"    Irrep    Number of SO\n");
  fprintf(outfile,"    -----    ------------\n");
  for(i=0;i<nirreps;i++)
    fprintf(outfile,"    %3d      %7d\n",i+1,num_so_per_irrep[i]);
  fprintf(outfile,"\n");
  if (print_lvl >= DEBUGPRINT) {
    fprintf(outfile,"    Prim#     Exponent     Norm. contr. coeff.\n");
    fprintf(outfile,"    -----  --------------  -------------------\n");
    for(i=0;i<num_prims;i++)
      fprintf(outfile,"    %3d    %14.7lf    %15.10lf\n",i+1,exponents[i],contr_coeff[i]);
    fprintf(outfile,"\n");
    fprintf(outfile,"    Shell#  Nuc#  L  SPRIM  SLOC  SNUMG\n");
    fprintf(outfile,"    ------  ----  -  -----  ----  -----\n");
    for(i=0;i<num_shells;i++)
      fprintf(outfile,"    %4d    %3d  %2d  %3d    %3d   %3d\n",i+1,
	      shell_nucleus[i]+1,shell_ang_mom[i],first_prim_shell[i]+1,first_basisfn_shell[i]+1,nprim_in_shell[i]);
    fprintf(outfile,"\n");
    fprintf(outfile,"    AMOrd. Shell#  Canon. Shell#\n");
    fprintf(outfile,"    -------------  -------------\n");
    for(i=0;i<num_shells;i++)
      fprintf(outfile,"    %7d        %7d\n",i+1,am2canon_shell_order[i]+1);
    fprintf(outfile,"\n\n");
  }
  fflush(outfile);

  return;
}


void cleanup()
{
  int i,l,uc,class,class_first,class_last;
  
  /*----------------------
    Free up global arrays
   ----------------------*/

/*  free_char_matrix(elem_name,NUM_ELEMENTS);*/
  free(sym_oper);
  free_block(geometry);
  free_matrix(full_geom,num_entries);
  free_block(Rref);
  free(nuclear_charges);
/*  free_char_matrix(element,num_atoms);
    free_char_matrix(full_element,num_entries);
  free_char_matrix(atom_basis,num_atoms);*/
  free_int_matrix(atom_orbit,num_atoms);
  free_int_matrix(class_orbit,num_classes);
  free_int_matrix(red_unique_orbit,num_uniques);
  for(uc=0;uc<num_unique_classes;uc++) {
    class = uc2c[uc];
    class_first = class;
    class_last = class + unique_class_degen[uc];
    for(i=class_first;i<class_last;i++)
      for(l=0;l<=max_angmom_class[i];l++) {
	free_matrix(class_so_coeff[i][l],ioff[l+1]*unique_class_degen[uc]);
	free(num_cart_so_in_class[i][l]);
	if (puream)
	  free(num_pureang_so_in_class[i][l]);
      }
  }
  for(i=0;i<num_classes;i++) {
    free(class_so_coeff[i]);
    free(num_cart_so_in_class[i]);
    if (puream)
      free(num_pureang_so_in_class[i]);
  }
  free(class_so_coeff);
  free(num_cart_so_in_class);
  if (puream)
    free(num_pureang_so_in_class);
  free(max_angmom_class);
  free(u2a);
  free(uc2c);
  free(unique_degen);
  free(unique_class_degen);
  free(num_cart_so_per_irrep);
  free_int_matrix(num_cart_so,MAX(max_angmom+1,MAXANGMOM));
  if (puream) {
    free(num_so_per_irrep);
    free_int_matrix(num_pureang_so,max_angmom+1);
    free_int_matrix(num_redun_so,max_angmom+1);
  }
  free(nshells_per_atom);
  free(first_shell_on_atom);
  free(exponents);
  free(contr_coeff);
  free(first_prim_shell);
  free(shell_nucleus);
  free(shell_ang_mom);
  free(nprim_in_shell);
  free(first_ao_shell);
  free(shells_per_am);
  free(am2canon_shell_order);
  for(i=0;i<MAX(max_angmom+1,MAXANGMOM);i++)
    free_matrix(ao_type_transmat[i],nirreps);
  free(ao_type_transmat);
  if (puream) {
    free(pureang_so_l);
    free(pureang_so_m);
    for(i=0;i<=max_angmom;i++)
      free_matrix(cart2pureang[i],2*i+1);
    free(cart2pureang);
  }
  free_matrix(usotao,num_so);
  free(atom_class);
  free(atom_position);
  free_int_matrix(ao_type_irr,MAX(max_angmom+1,MAXANGMOM));
  free(ioff);
  free(df);
  free(z_geom);
  
  return;
}
