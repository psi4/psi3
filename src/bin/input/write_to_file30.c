#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <math.h>
#include <file30_params.h>
#include "input.h"
#include "global.h"
#include "defines.h"


/*-----------------------------------------------------------------------------------------------------------------
  This function writes information out to file30
 -----------------------------------------------------------------------------------------------------------------*/

void write_to_file30(double repulsion)
{

  PSI_FPTR ptr, junk;
  int *constants,*pointers,*calcs;
  char *calc_label;
  int i,j,k,l;
  int atom,class,ua,shell,irr,symop,shell_first,shell_last,us;
  int eq_atom;
  int *arr_int;
  double *arr_double;
  int *ict;
  double *cspd;     /*Array of contraction coefficients in file30 format*/
  char *atom_label;
  int **shell_transm;
  int **ioff_irr;
  int max_num_prims;
  int max_atom_degen;
  int max_angmom_unique;
  
  
  constants = init_int_array(MCONST);
  pointers = init_int_array(MPOINT);
  calcs = init_int_array(MPOINT);

  /* Check the max_angmom. If it's > 4 - die */
  if (max_angmom >= MAXANGMOM)
    punt("Angular momentum is too high to be handled by file30");

  
  /*----------------------------------
    Write out the label then 80 zeros
   ----------------------------------*/
  rfile(CHECKPOINTFILE);
  calc_label = init_char_array(80);
  strncpy(calc_label,label,MIN(80,strlen(label)));
  wwritw(CHECKPOINTFILE,(char *) calc_label, 20*(sizeof(int)),0, &ptr);
  wwritw(CHECKPOINTFILE,(char *) pointers, 80*(sizeof(int)),ptr,&ptr);


  /*------------------------------------------------------
    Zero arrays for constants, pointers, and calculations
   ------------------------------------------------------*/
  wwritw(CHECKPOINTFILE,(char *) constants, MCONST*(sizeof(int)),100*sizeof(int),&ptr);
  wwritw(CHECKPOINTFILE,(char *) pointers, MPOINT*(sizeof(int)),300*sizeof(int),&ptr);
  wwritw(CHECKPOINTFILE,(char *) calcs, MCALCS*(sizeof(int)),500*sizeof(int),&ptr);


  /*-----------------------------------
    Start writing pointers to the file
   -----------------------------------*/

  /* Nuclear charges */
  pointers[0] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) nuclear_charges, num_atoms*(sizeof(double)),ptr,&ptr);

  /* Transformation table for atoms - just atom_orbit transposed */
  pointers[1] = ptr/sizeof(int) + 1;
  ict = init_int_array(num_atoms);
  for(i=0;i<nirreps;i++) {
    for(j=0;j<num_atoms;j++)
      ict[j] = atom_orbit[j][i]+1;
    wwritw(CHECKPOINTFILE,(char *) ict, num_atoms*(sizeof(int)),ptr,&ptr);
  }
  free(ict);

  /* Number of shells per atom */
  pointers[2] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) nshells_per_atom, num_atoms*(sizeof(int)),ptr,&ptr);

  /* Pointer to shells on atoms */
  pointers[3] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_atoms);
  for(i=0;i<num_atoms;i++)
    arr_int[i] = first_shell_on_atom[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_atoms*(sizeof(int)),ptr,&ptr);
  free(arr_int);

  /* Exponents of primitive gaussians */
  pointers[4] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) exponents, num_prims*(sizeof(double)),ptr,&ptr);

  /* Contraction coefficients */
  pointers[5] = ptr/sizeof(int) + 1;
     /*------This piece of code is for segmented contractions ONLY------*/
  cspd = init_array(MRCRU*num_prims*MAXANGMOM);
  for(j=0;j<num_shells;j++)
    for(k=0;k<nprim_in_shell[j];k++)
	/*---
	  Pitzer normalization of Psi 2 is NOT used - cc's for d-functions used to be
	  multiplied by sqrt(3), f - by sqrt(15), g - sqrt(105), etc
	  ---*/
      cspd[shell_ang_mom[j]*num_prims+first_prim_shell[j]+k] = contr_coeff[first_prim_shell[j]+k];
  wwritw(CHECKPOINTFILE,(char *) cspd, MRCRU*num_prims*MAXANGMOM*(sizeof(double)),ptr,&ptr);
  free(cspd);

  /* Pointer to primitives for a shell */
  pointers[6] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_prim_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);

  /* Atom on which nth shell is centered */
  pointers[7] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = shell_nucleus[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);

  /* Angular momentum of a shell */
  pointers[8] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = shell_ang_mom[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);

  /* Number of contracted functions (primitives) in a shell */
  pointers[9] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) nprim_in_shell, num_shells*(sizeof(int)),ptr,&ptr);

  /* Pointer to the first AO in shell */
  pointers[10] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_ao_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  free(arr_int);

  /* Labels of irreps */
  pointers[15] = ptr/sizeof(int) + 1;
  for(i=0;i<nirreps;i++)
    wwritw(CHECKPOINTFILE,(char *) irr_labels[i], 4,ptr,&ptr);

  /* Class an atom belongs to*/
  pointers[20] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_atoms);
  for(atom=0;atom<num_atoms;atom++)
    arr_int[atom] = atom_class[atom] + 1;
/*  k = arr_int[1]; arr_int[1] = arr_int[2]; arr_int[2] = k;*/
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_atoms*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* Pointer to start of nth symmetry block*/
  pointers[22] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(nirreps);
  arr_int[0] = 0;
  for(irr=1;irr<nirreps;irr++)
    arr_int[irr] = arr_int[irr-1] + ioff[num_cart_so_per_irrep[irr-1]];
  wwritw(CHECKPOINTFILE,(char *) arr_int, nirreps*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* Number of cartesian SO's of nth symmetry */
  pointers[23] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) num_cart_so_per_irrep, nirreps*sizeof(int),ptr,&ptr);

  /* Transformation matrices for coordinates (or p-functions) */
  pointers[24] = ptr/sizeof(int) + 1;
  arr_double = init_array(3);
  for(symop=0;symop<nirreps;symop++)
    for(i=0;i<3;i++) {
      arr_double[i] = ao_type_transmat[1][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 3*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }
  free(arr_double);

  /* Transformation matrix for shells */
  pointers[26] = ptr/sizeof(int) + 1;
  shell_transm = init_int_matrix(num_shells,nirreps);;
  for(atom=0;atom<num_atoms;atom++)
    for(symop=0;symop<nirreps;symop++) {
      eq_atom = atom_orbit[atom][symop];
      for(i=0;i<nshells_per_atom[atom];i++)
	shell_transm[first_shell_on_atom[atom]+i][symop] = first_shell_on_atom[eq_atom] + i + 1;
    }
  for(i=0;i<num_shells;i++)
    wwritw(CHECKPOINTFILE,(char *) shell_transm[i], nirreps*(sizeof(int)),ptr,&ptr);
  free_int_matrix(shell_transm,num_shells);

  /* Labels of atoms */
  pointers[27] = ptr/sizeof(int) + 1;
  atom_label = init_char_array(8);
  for(atom=0;atom<num_atoms;atom++) {
    strncpy(atom_label,element[atom],strlen(element[atom]));
    wwritw(CHECKPOINTFILE,(char *) atom_label, 8*(sizeof(char)),ptr,&ptr);
  }
  
  free(atom_label);

  /* Orbitals per irrep */
  pointers[36] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) num_so_per_irrep, nirreps*sizeof(int),ptr,&ptr);  

  /* Symmetry label */
  pointers[37] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) symmetry, 4*sizeof(char),ptr,&ptr);

  /* Symmetry positions of atoms - for more info see count_uniques.c */
  pointers[38] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) atom_position, num_atoms*sizeof(int),ptr,&ptr);

  /* Unique shell number to full shell number mapping array */
  pointers[39] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_unique_shells);
  us = 0;
  for(ua=0;ua<num_uniques;ua++) {
    atom = u2a[ua];
    shell = first_shell_on_atom[atom];
    for(i=0;i<nshells_per_atom[atom];i++,shell++,us++)
      arr_int[us] = shell;
  }
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_unique_shells*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* SO to AO transformation matrix */
  pointers[40] = ptr/sizeof(int) + 1;
  for(i=0;i<num_so;i++)
    wwritw(CHECKPOINTFILE,(char *) usotao[i], num_ao*sizeof(double),ptr,&ptr);

  /* SO to basis functions transformation matrix */
  if (puream) {
  pointers[41] = ptr/sizeof(int) + 1;
  for(i=0;i<num_so;i++)
    wwritw(CHECKPOINTFILE,(char *) usotbf[i], num_so*sizeof(double),ptr,&ptr);
  }

  /* Pointers to first basis functions from shells */
  pointers[42] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_basisfn_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  free(arr_int);

  /* Unique atom number to full atom number mapping array */
  pointers[43] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) u2a, num_uniques*sizeof(int),ptr,&ptr);

  /* Mapping between canonical Cotton ordering of symmetry operations
     in the point group to the symmetry.h-defined ordering */
  pointers[44] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) sym_oper, nirreps*sizeof(int),ptr,&ptr);

  
  /*---------------------------
    Write pointers to the file
   ---------------------------*/

  wwritw(CHECKPOINTFILE,(char *) pointers, MPOINT*sizeof(int),300*sizeof(int),&junk);


  /*-------------------------------------
    Write directory to calcs to the file
   -------------------------------------*/

  calcs[0] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) calcs, MCALCS*(sizeof(int)),500*sizeof(int),&junk);
  wwritw(CHECKPOINTFILE,(char *) label, 20*(sizeof(int)),ptr, &ptr);
  arr_int = init_int_array(40);
  wwritw(CHECKPOINTFILE,(char *) arr_int, 40*sizeof(int),ptr,&ptr);
  /* Zero pointers to vectors */
  wwritw(CHECKPOINTFILE,(char *) arr_int, 20*sizeof(int),ptr,&ptr);
  free(arr_int);
  /* Write out geometry */
  for(atom=0;atom<num_atoms;atom++)
    wwritw(CHECKPOINTFILE,(char *) geometry[atom], 3*sizeof(double),ptr,&ptr);
  /* Zero energies */
  arr_double = init_array(10);
  arr_double[0] = repulsion;
  wwritw(CHECKPOINTFILE,(char *) arr_double, 10*sizeof(double),ptr,&ptr);
  free(arr_double);
  fprintf(outfile,"    Wrote %u bytes to FILE%d\n\n",ptr,CHECKPOINTFILE);
  
  /*-----------------------------------
    Put all constants into constants[]
   -----------------------------------*/

  constants[0] = ptr/sizeof(int) + 1;
  constants[1] = MPOINT;
  constants[2] = MCONST;
  constants[3] = MCALCS;
  constants[4] = 1;
  constants[5] = num_uniques;
  constants[6] = (int) rotor;
  constants[7] = num_unique_shells;
  constants[17] = num_so;
  constants[18] = num_atoms;
  constants[21] = num_ao;
  constants[26] = num_shells;
  constants[27] = nirreps;
  constants[31] = num_prims;
  /*--- These are SCF constants to be filled by CSCF ---*/
  constants[40] = 0;
  constants[41] = 0;
  constants[42] = 0;
  constants[44] = 0;
  constants[45] = 0;
  constants[46] = 0;
  constants[50] = 0;
  wwritw(CHECKPOINTFILE,(char *) constants, MCONST*sizeof(int),100*sizeof(int),&junk);

  rclose(CHECKPOINTFILE,3);
  free(constants);
  free(pointers);
  free(calcs);
  return;
}
