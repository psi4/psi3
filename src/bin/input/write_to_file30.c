#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <math.h>
#include <file30_params.h>
#include "input.h"
#include "global.h"
#include "defines.h"


void pack_4int(int **, int *, int, int);

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
  int **isc, **ipc, **iso;
  int **loc, **loc2;
  int **ioff_irr;
  int max_num_prims;
  int max_atom_degen;
  int max_angmom_unique;
  int mru,kaords;
  
  
  constants = init_int_array(MCONST);
  pointers = init_int_array(MPOINT);
  calcs = init_int_array(MPOINT);

  /* Check the max_angmom. If it's > 4 - die */
  if (max_angmom >= MAXANGMOM)
    punt("  You went too high on the angular momentum stair!\n\n");

  
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
	/*-------Pitzer normalization is used - cc's for d-functions are multiplied by sqrt(3), f - by sqrt(15), g - sqrt(105), etc */
      cspd[shell_ang_mom[j]*num_prims+first_prim_shell[j]+k] = contr_coeff[first_prim_shell[j]+k]/**sqrt(df[2*shell_ang_mom[j]])*/;
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

  /* First ao type in shell (S,X,Y,Z,XX,..) */
  pointers[11] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_ao_type_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);

  /* Last ao type in shell (S,X,Y,Z,XX,..) */
  pointers[12] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = last_ao_type_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  free(arr_int);

  /* Packed SO transformation matrices - for abelian point group it's just an array of 1's*/
  pointers[13] = ptr/sizeof(int) + 1;
  isc = init_int_matrix(num_shells,nirreps);
  for(i=0;i<num_shells;i++)
    for(j=0;j<nirreps;j++)
      isc[i][j] = 1;
  arr_int = (int *) malloc(num_shells*((int)(nirreps+3)/4)*sizeof(int));
  pack_4int(isc,arr_int,num_shells,nirreps);
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*((int)(nirreps+3)/4)*(sizeof(int)),ptr,&ptr);

  /* Packed SO transformation information */
  pointers[14] = ptr/sizeof(int) + 1;
  ipc = isc;
  ioff_irr = init_int_matrix(max_angmom+1,nirreps);
  for(class=0;class<num_classes;class++) {
    for(l=0;l<=max_angmom;l++)
      for(irr=0;irr<nirreps;irr++)
	ioff_irr[l][irr] = 1;
    for(atom=0;atom<num_atoms;atom++)
      if (atom_class[atom] == class) {
	shell_first = first_shell_on_atom[atom];
	shell_last = first_shell_on_atom[atom] + nshells_per_atom[atom];
	for(shell=shell_first;shell<shell_last;shell++) {
	  l = shell_ang_mom[shell];
	  for(irr=0;irr<nirreps;irr++) {
	    ipc[shell][irr] = ioff_irr[l][irr];
	    ioff_irr[l][irr] += num_cart_so_in_class[class][l][irr];
	  }
	}
      }
  }
  free_int_matrix(ioff_irr,max_angmom+1);
  pack_4int(ipc,arr_int,num_shells,nirreps);
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*((int)(nirreps+3)/4)*(sizeof(int)),ptr,&ptr);
  free(arr_int);

  /* Labels of irreps */
  pointers[15] = ptr/sizeof(int) + 1;
  for(i=0;i<nirreps;i++)
    wwritw(CHECKPOINTFILE,(char *) irr_labels[i], 4,ptr,&ptr);

  /* Coefficients of SO's */
  pointers[16] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) angso_coeff, num_angso_coeff*sizeof(double),ptr,&ptr);

  /* Labels of SO matrices */
  pointers[17] = ptr/sizeof(int) + 1;
  arr_int = (int *) malloc(num_angso_coeff*sizeof(int));
  pack_4int(angso_labels,arr_int,num_angso_coeff,4);
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_angso_coeff*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* Pointers to beginning of SO matrix */
  pointers[18] = ptr/sizeof(int) + 1;
  loc = init_int_matrix(MAXANGMOM,num_atoms);
  loc2 = init_int_matrix(MAXANGMOM,num_atoms);
  for(l=0;l<MAXANGMOM;l++)
    for(atom=0;atom<num_atoms;atom++)
      loc[l][atom] = -1;
  for(class=0;class<num_classes;class++)
    for(atom=0;atom<num_atoms;atom++)
      if (atom_class[atom] == class) {
	for(l=0;l<=max_angmom_class[class];l++) {
	  loc[l][atom] = first_angso_coeff_class[class][l];
	  loc2[l][atom] = last_angso_coeff_class[class][l];
	}
	break;
      }
  for(l=0;l<MAXANGMOM;l++)
    wwritw(CHECKPOINTFILE,(char *) loc[l], num_atoms*sizeof(int),ptr,&ptr);
  free_int_matrix(loc,MAXANGMOM);

  /* Pointers to end of SO matrix */
  pointers[19] = ptr/sizeof(int) + 1;
  for(l=0;l<MAXANGMOM;l++)
    wwritw(CHECKPOINTFILE,(char *) loc2[l], num_atoms*sizeof(int),ptr,&ptr);
  free_int_matrix(loc2,MAXANGMOM);

  /* Class an atom belongs to*/
  pointers[20] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_atoms);
  for(atom=0;atom<num_atoms;atom++)
    arr_int[atom] = atom_class[atom] + 1;
/*  k = arr_int[1]; arr_int[1] = arr_int[2]; arr_int[2] = k;*/
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_atoms*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* Degeneracy of nth irrep*/
  pointers[21] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(nirreps);
  for(irr=0;irr<nirreps;irr++)
    arr_int[irr] = 1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, nirreps*sizeof(int),ptr,&ptr);

  /* Pointer to start of nth symmetry block*/
  pointers[22] = ptr/sizeof(int) + 1;
  arr_int[0] = 0;
  for(irr=1;irr<nirreps;irr++)
    arr_int[irr] = arr_int[irr-1] + ioff[num_cart_so_per_irrep[irr-1]];
  wwritw(CHECKPOINTFILE,(char *) arr_int, nirreps*sizeof(int),ptr,&ptr);

  /* Number of cartesian SO's of nth symmetry */
  pointers[23] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) num_cart_so_per_irrep, nirreps*sizeof(int),ptr,&ptr);

  /* Transformation matrices for coordinates */
  pointers[24] = ptr/sizeof(int) + 1;
  arr_double = init_array(3);
  for(symop=0;symop<nirreps;symop++)
    for(i=0;i<3;i++) {
      arr_double[i] = ao_type_transmat[1][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 3*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }
  free(arr_double);

  /* Inverse of symmetry operation */
  pointers[25] = ptr/sizeof(int) + 1;
  for(symop=0;symop<nirreps;symop++)
    arr_int[symop] = symop + 1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, nirreps*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* Packed transformation matrix for shells */
  pointers[26] = ptr/sizeof(int) + 1;
  iso = isc;
  for(atom=0;atom<num_atoms;atom++)
    for(symop=0;symop<nirreps;symop++) {
      eq_atom = atom_orbit[atom][symop];
      for(i=0;i<nshells_per_atom[atom];i++)
	iso[first_shell_on_atom[atom]+i][symop] = first_shell_on_atom[eq_atom] + i + 1;
    }
  arr_int = (int *) malloc(num_shells*((int)(nirreps+3)/4)*sizeof(int));
  pack_4int(iso,arr_int,num_shells,nirreps);
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*((int)(nirreps+3)/4)*(sizeof(int)),ptr,&ptr);
  free(arr_int);
  free_int_matrix(isc,num_shells);

  /* Labels of atoms */
  pointers[27] = ptr/sizeof(int) + 1;
  atom_label = init_char_array(8);
  for(atom=0;atom<num_atoms;atom++) {
    strncpy(atom_label,element[atom],strlen(element[atom]));
    wwritw(CHECKPOINTFILE,(char *) atom_label, 8*(sizeof(char)),ptr,&ptr);
  }
  
  free(atom_label);

  /* Transformation matrices for p-functions */
  pointers[30] = ptr/sizeof(int) + 1;
  arr_double = init_array(ioff[MAXANGMOM+1]);
  for(symop=0;symop<nirreps;symop++)
    for(i=0;i<3;i++) {
      arr_double[i] = ao_type_transmat[1][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 3*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }

  /* Transformation matrices for d-functions */
  pointers[31] = ptr/sizeof(int) + 1;
  for(symop=0;symop<nirreps;symop++)
    for(i=0;i<6;i++) {
      arr_double[i] = ao_type_transmat[2][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 6*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }

  /* Transformation matrices for f-functions */
  pointers[32] = ptr/sizeof(int) + 1;
  for(symop=0;symop<nirreps;symop++)
    for(i=0;i<10;i++) {
      arr_double[i] = ao_type_transmat[3][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 10*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }

  /* Transformation matrices for g-functions */
  pointers[33] = ptr/sizeof(int) + 1;
  for(symop=0;symop<nirreps;symop++)
    for(i=0;i<15;i++) {
      arr_double[i] = ao_type_transmat[4][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 15*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }
  free(arr_double);

  /* Number of contractions in each shell - 1 if only segmented basis sets are used */
  pointers[34] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = 1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*sizeof(int),ptr,&ptr);
  free(arr_int);

  /* Shell type - 1 if pure angular momentum, 0 if cartesian */
  pointers[35] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(MAXANGMOM);
  arr_int[0] = 0;
  if (max_angmom >= 1)
    arr_int[1] = 0;
  if (puream)
    for(l=2;l<=max_angmom;l++)
      arr_int[l] = 1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, MAXANGMOM*sizeof(int),ptr,&ptr);
  free(arr_int);

  /*----------------------------------------------
    Experimental additions - needed to feed CINTS
   ----------------------------------------------*/
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

  /* Irrep character rows packed in int's */
  pointers[43] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) irr_char_str, nirreps*sizeof(int),ptr,&ptr);

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
  /* Compute some weird constants */
  max_num_prims = 0;
  max_atom_degen = 1;
  mru = 0;
  kaords = 0;
  for(ua=0;ua<num_uniques;ua++) {
    max_atom_degen = MAX(max_atom_degen,unique_degen[ua]);
    atom = u2a[ua];
    shell_first = first_shell_on_atom[atom];
    shell_last = first_shell_on_atom[atom] + nshells_per_atom[atom];
    max_angmom_unique = 0; 
    for(shell=shell_first;shell<shell_last;shell++) {
      max_angmom_unique = MAX(max_angmom_unique,shell_ang_mom[shell]);
      max_num_prims = MAX(max_num_prims,nprim_in_shell[shell]);
    }
    kaords += max_angmom_unique + 1;
    mru = MAX(mru,unique_degen[ua]*ioff[max_angmom_unique+1]);
  }
      
  constants[6] = kaords;
  constants[7] = num_unique_shells;
  constants[8] = kaords;
  constants[9] = nirreps;
  constants[10] = mru;
  constants[11] = mru;
  constants[12] = mru;
  constants[13] = max_num_prims;
  constants[14] = max_atom_degen;
  constants[15] = num_unique_shells;
  constants[16] = nirreps;
  constants[17] = num_so;
  constants[18] = num_atoms;
  constants[19] = 0;
  constants[20] = 0;
  constants[21] = num_ao;
  constants[22] = 0;
  constants[23] = num_so*(num_so+1)/2;
  constants[24] = 0;
  constants[25] = 0;
  constants[26] = num_shells;
  constants[27] = nirreps;
  constants[28] = nirreps;
  constants[29] = 40;
  constants[30] = 48;
  constants[31] = num_prims;
  constants[32] = 100;
  constants[33] = 14;
  constants[34] = num_angso_coeff;
  constants[35] = 9*nirreps;
  constants[36] = 12;
  constants[37] = 4;
  constants[38] = (nirreps+3)/4;
  constants[39] = (nirreps+3)/4;
  constants[40] = 0;
  constants[41] = 0;
  constants[42] = 0;
  constants[43] = MRCRU;
  wwritw(CHECKPOINTFILE,(char *) constants, MCONST*sizeof(int),100*sizeof(int),&junk);

  rclose(CHECKPOINTFILE,3);
  free(constants);
  free(pointers);
  free(calcs);
  return;
}



/* pack_4int()   Packs a matrix of small integers (normally, symmetry information)
**               into a sequential array of integers.
**
**  FORMAT: source is converted into a packed form. Each element of a row is
**  represented by a byte, therefore if nirreps <= 4 - each row is packed   
**  into an integer word, else if nirreps == 8 - into two integer words.
**   Here's how it works :
**   | source[3] | source[2] | source[1] | source[0] |                 
**   leftmost byte                       rightmost byte
**			
**
**  arguments: int **source - integer matrix num_rows by nirreps
**             int  *dest   - the destination array of
**                            num_rows*((nirreps == 8) ? 2 : 1) integers
**             int num_rows - number of rows in source
**             int nirreps  - number of columns in source
**
**  returns: nothing	
*/

void pack_4int(int **source, int *dest, int num_rows, int num_cols)
{
  int i;

  switch (num_cols) {
    case 8: bzero((char *) (dest+num_rows),sizeof(int)*num_rows);
    case 4:
    case 2:
    case 1: bzero((char *) dest,sizeof(int)*num_rows);
  }
  for(i=0;i<num_rows;i++)
    switch (num_cols) {
      case 8: dest[i+num_rows] += source[i][4];
              dest[i+num_rows] += (source[i][5] << 8);
              dest[i+num_rows] += (source[i][6] << 16);
              dest[i+num_rows] += (source[i][7] << 24);
      case 4: dest[i] += (source[i][2] << 16);
              dest[i] += (source[i][3] << 24);
      case 2: dest[i] += (source[i][1] << 8);
      case 1: dest[i] += source[i][0];
    }

}
