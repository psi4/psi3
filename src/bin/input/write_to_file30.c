#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <file30_params.h>
#include <psifiles.h>
#include "input.h"
#include "global.h"
#include "defines.h"

void write_scf_calc(PSI_FPTR *ptr, int *pointers);

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
  double **mat_double;
  int *ict, **ict_tmp;
  int *scf_pointers;
  double *cspd;     /*Array of contraction coefficients in file30 format*/
  char *atom_label, **tmp_atom_label;
  int **shell_transm;
  int **ioff_irr;
  int max_num_prims;
  int max_atom_degen;
  int max_angmom_unique;
  psio_address chkptr;
  
  constants = init_int_array(MCONST);
  pointers = init_int_array(MPOINT);
  calcs = init_int_array(MPOINT);

  /* Check the max_angmom. If it's >= MAXANGMOM from file30_params.h - die */
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

  /*----------------------------------
    Write out basic info to chkpt
    ---------------------------------*/
  chkpt_init();
  chkpt_wt_label(calc_label);
  chkpt_wt_num_unique_atom(num_uniques);
  chkpt_wt_num_unique_shell(num_unique_shells);
  chkpt_wt_rottype(rotor);
  chkpt_wt_max_am(max_angmom);
  chkpt_wt_nso(num_so);
  chkpt_wt_nao(num_ao);
  chkpt_wt_nshell(num_shells);
  chkpt_wt_nirreps(nirreps);
  chkpt_wt_nprim(num_prims);
  chkpt_wt_natom(num_atoms);
  chkpt_wt_nentry(num_entries);

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
  chkpt_wt_zvals(nuclear_charges);

  /* Transformation table for atoms - just atom_orbit transposed */
  pointers[1] = ptr/sizeof(int) + 1;
  ict = init_int_array(num_atoms);
  ict_tmp = init_int_matrix(nirreps, num_atoms);
  for(i=0;i<nirreps;i++) {
    for(j=0;j<num_atoms;j++) {
      ict[j] = atom_orbit[j][i]+1;
      ict_tmp[i][j] = ict[j];
    }
    wwritw(CHECKPOINTFILE,(char *) ict, num_atoms*(sizeof(int)),ptr,&ptr);
  }
  chkpt_wt_ict(ict_tmp);
  free_int_matrix(ict_tmp, nirreps);
  free(ict);

  /** NOT ACTUALLY USED **/
  /* Number of shells per atom */
  pointers[2] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) nshells_per_atom, num_atoms*(sizeof(int)),ptr,&ptr);

  /** NOT ACTUALLY USED **/
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
  chkpt_wt_exps(exponents);

  /* Contraction coefficients */
  pointers[5] = ptr/sizeof(int) + 1;
  /*------This piece of code is for segmented contractions ONLY------*/
  cspd = init_array(num_prims*MAXANGMOM);
  for(j=0;j<num_shells;j++)
    for(k=0;k<nprim_in_shell[j];k++)
      /*---
	Pitzer normalization of Psi 2 is NOT used - cc's for d-functions used to be
	multiplied by sqrt(3), f - by sqrt(15), g - sqrt(105), etc
	---*/
      cspd[shell_ang_mom[j]*num_prims+first_prim_shell[j]+k] = contr_coeff[first_prim_shell[j]+k];
  wwritw(CHECKPOINTFILE,(char *) cspd, num_prims*MAXANGMOM*(sizeof(double)),ptr,&ptr);
  chkpt_wt_contr(cspd);
  free(cspd);

  /* Pointer to primitives for a shell */
  pointers[6] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_prim_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  chkpt_wt_sprim(arr_int);

  /* Atom on which nth shell is centered */
  pointers[7] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = shell_nucleus[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  chkpt_wt_snuc(arr_int);

  /* Angular momentum of a shell */
  pointers[8] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = shell_ang_mom[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  chkpt_wt_stype(arr_int);

  /* Number of contracted functions (primitives) in a shell */
  pointers[9] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) nprim_in_shell, num_shells*(sizeof(int)),ptr,&ptr);
  chkpt_wt_snumg(nprim_in_shell);

  /* Pointer to the first AO in shell */
  pointers[10] = ptr/sizeof(int) + 1;
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_ao_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  chkpt_wt_sloc(arr_int);
  free(arr_int);

  /* Labels of irreps */
  pointers[15] = ptr/sizeof(int) + 1;
  for(i=0;i<nirreps;i++) {
    wwritw(CHECKPOINTFILE,(char *) irr_labels[i], 4,ptr,&ptr);
  }
  chkpt_wt_irr_labs(irr_labels);

  /*** NOT APPARENTLY USED ***/
  /* Class an atom belongs to*/
  pointers[20] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_atoms);
  for(atom=0;atom<num_atoms;atom++)
    arr_int[atom] = atom_class[atom] + 1;
  /*  k = arr_int[1]; arr_int[1] = arr_int[2]; arr_int[2] = k;*/
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_atoms*sizeof(int),ptr,&ptr);
  free(arr_int);

  /*** NOT APPARENTLY USED ***/
  /* Pointer to start of nth symmetry block*/
  pointers[22] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(nirreps);
  arr_int[0] = 0;
  for(irr=1;irr<nirreps;irr++)
    arr_int[irr] = arr_int[irr-1] + ioff[num_cart_so_per_irrep[irr-1]];
  wwritw(CHECKPOINTFILE,(char *) arr_int, nirreps*sizeof(int),ptr,&ptr);
  free(arr_int);

  /*** NOT APPARENTLY USED ***/
  /* Number of cartesian SO's of nth symmetry */
  pointers[23] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) num_cart_so_per_irrep, nirreps*sizeof(int),ptr,&ptr);

  /*** NOT APPARENTLY USED ***/
  /* Transformation matrices for coordinates (or p-functions) */
  pointers[24] = ptr/sizeof(int) + 1;
  arr_double = init_array(3);
  mat_double = block_matrix(nirreps,9);
  for(symop=0;symop<nirreps;symop++) {
    for(i=0;i<3;i++) {
      arr_double[i] = ao_type_transmat[1][symop][i];
      wwritw(CHECKPOINTFILE,(char *) arr_double, 3*sizeof(double),ptr,&ptr);
      arr_double[i] = 0.0;
    }
    mat_double[symop][0] = ao_type_transmat[1][symop][0];
    mat_double[symop][4] = ao_type_transmat[1][symop][1];
    mat_double[symop][8] = ao_type_transmat[1][symop][2];
  }
  free(arr_double);
  chkpt_wt_cartrep(mat_double);
  free(mat_double);

  /* Transformation matrix for shells */
  pointers[26] = ptr/sizeof(int) + 1;
  shell_transm = init_int_matrix(num_shells,nirreps);;
  for(atom=0;atom<num_atoms;atom++)
    for(symop=0;symop<nirreps;symop++) {
      eq_atom = atom_orbit[atom][symop];
      for(i=0;i<nshells_per_atom[atom];i++)
	shell_transm[first_shell_on_atom[atom]+i][symop] = first_shell_on_atom[eq_atom] + i + 1;
    }
  for(i=0;i<num_shells;i++) {
    wwritw(CHECKPOINTFILE,(char *) shell_transm[i], nirreps*(sizeof(int)),ptr,&ptr);
  }
  chkpt_wt_shell_transm(shell_transm);
  free_int_matrix(shell_transm,num_shells);

  /*** NOT APPARENTLY USED ***/
  /* Labels of atoms */
  pointers[27] = ptr/sizeof(int) + 1;
  atom_label = init_char_array(8);
  for(atom=0;atom<num_atoms;atom++) {
    strncpy(atom_label,element[atom],strlen(element[atom]));
    wwritw(CHECKPOINTFILE,(char *) atom_label, 8*(sizeof(char)),ptr,&ptr);
  }
  free(atom_label);

  /* Labels of atoms including dummy atoms */
  if(cartOn)
    j = num_atoms;
  else j = num_entries;
  pointers[28] = ptr/sizeof(int) + 1;
  atom_label = init_char_array(8);
  for(atom=0;atom<j;atom++) {
    strncpy(atom_label,full_element[atom],strlen(full_element[atom]));
    atom_label[strlen(full_element[atom])] = '\0';
    wwritw(CHECKPOINTFILE,(char *) atom_label, 8*(sizeof(char)),ptr,&ptr);
  }
  free(atom_label);

  tmp_atom_label = (char **) malloc(j*sizeof(char *));
  for(atom=0; atom<j; atom++) {
    tmp_atom_label[atom] = init_char_array(MAX_ELEMNAME);

    if(strlen(full_element[atom]) > MAX_ELEMNAME)
      punt("Element name exceeds limit, MAX_ELEMNAME.\n");

    strncpy(tmp_atom_label[atom],full_element[atom],strlen(full_element[atom]));
    tmp_atom_label[atom][strlen(full_element[atom])] = '\0';
  }
  chkpt_wt_felement(tmp_atom_label);
  for(atom=0; atom<j; atom++) {
    free(tmp_atom_label[atom]);
  }
  free(tmp_atom_label);

  /* Orbitals per irrep */
  pointers[36] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) num_so_per_irrep, nirreps*sizeof(int),ptr,&ptr);  
  chkpt_wt_sopi(num_so_per_irrep);

  /* Symmetry label */
  pointers[37] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) symmetry, 4*sizeof(char),ptr,&ptr);
  chkpt_wt_sym_label(symmetry);

  /* Symmetry positions of atoms - for more info see count_uniques.c */
  pointers[38] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) atom_position, num_atoms*sizeof(int),ptr,&ptr);
  chkpt_wt_atom_position(atom_position);

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
  chkpt_wt_us2s(arr_int);
  free(arr_int);

  /* SO to AO transformation matrix */
  pointers[40] = ptr/sizeof(int) + 1;
  for(i=0;i<num_so;i++) {
    wwritw(CHECKPOINTFILE,(char *) usotao[i], num_ao*sizeof(double),ptr,&ptr);
  }
  chkpt_wt_usotao(usotao);

  /* SO to basis functions transformation matrix */
  if (puream) {
    pointers[41] = ptr/sizeof(int) + 1;
    for(i=0;i<num_so;i++)
      wwritw(CHECKPOINTFILE,(char *) usotbf[i], num_so*sizeof(double),ptr,&ptr);
    chkpt_wt_usotbf(usotbf);
  }

  /* Pointers to first basis functions from shells */
  pointers[42] = ptr/sizeof(int) + 1;
  arr_int = init_int_array(num_shells);
  for(i=0;i<num_shells;i++)
    arr_int[i] = first_basisfn_shell[i]+1;
  wwritw(CHECKPOINTFILE,(char *) arr_int, num_shells*(sizeof(int)),ptr,&ptr);
  chkpt_wt_sloc_new(arr_int);

  free(arr_int);

  /* Unique atom number to full atom number mapping array */
  pointers[43] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) u2a, num_uniques*sizeof(int),ptr,&ptr);
  chkpt_wt_ua2a(u2a);

  /* Mapping between canonical Cotton ordering of symmetry operations
     in the point group to the symmetry.h-defined ordering */
  pointers[44] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) sym_oper, nirreps*sizeof(int),ptr,&ptr);
  chkpt_wt_symoper(sym_oper);

  /* write z_mat if it exists, see global.h for info about z_entry structure */
  if(!cartOn) {
    pointers[46] = ptr/sizeof(int) + 1;
    wwritw(CHECKPOINTFILE,(char *) z_geom, num_entries*(sizeof(struct z_entry)),ptr,&ptr);
    chkpt_wt_zmat(z_geom);
  }

  /* Number of shells in each angmom block */
  pointers[47] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) shells_per_am, (max_angmom+1)*sizeof(int),ptr,&ptr);
  chkpt_wt_shells_per_am(shells_per_am);

  /* Mapping array from the am-blocked to the canonical (in the order of
     appearance) ordering of shells */
  pointers[48] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) am2canon_shell_order, num_shells*sizeof(int),ptr,&ptr);
  chkpt_wt_am2canon_shell_order(am2canon_shell_order);

  /* Matrix representation of rotation back to the reference frame */
  pointers[49] = ptr/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) Rref[0], 9*sizeof(double),ptr,&ptr);
  chkpt_wt_rref(Rref);

  /* write full_geom, cartesian geometry with dummy atoms included */
  pointers[50] = ptr/sizeof(int) +1;
  for(i=0;i<num_entries;++i)
    wwritw(CHECKPOINTFILE,(char *) full_geom[i], 3*sizeof(double),ptr,&ptr);
  if(full_geom != NULL) chkpt_wt_fgeom(full_geom); /* only available if zmat input */

  
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
  free(arr_int);

  /* Zero pointers to vectors */
  scf_pointers = init_int_array(20);
  wwritw(CHECKPOINTFILE,(char *) scf_pointers, 20*sizeof(int),ptr,&ptr);
  free(scf_pointers);

  /* Write out geometry */
  for(atom=0;atom<num_atoms;atom++)
    wwritw(CHECKPOINTFILE,(char *) geometry[atom], 3*sizeof(double),ptr,&ptr);
  chkpt_wt_geom(geometry);

  /* Write out energies */
  arr_double = init_array(5);
  arr_double[0] = repulsion;
  wwritw(CHECKPOINTFILE,(char *) arr_double, 5*sizeof(double),ptr,&ptr);
  chkpt_wt_enuc(repulsion);
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
  constants[8] = max_angmom;
  constants[17] = num_so;
  constants[18] = num_atoms;
  if(cartOn)
    constants[19] = num_atoms;
  else
    constants[19] = num_entries;
  constants[21] = num_ao;
  constants[26] = num_shells;
  constants[27] = nirreps;
  constants[31] = num_prims;
  /*--- These are SCF constants ---*/
  constants[40] = 0;
  constants[41] = 0;
  constants[42] = 0;
  constants[44] = 0;
  constants[45] = 0;
  constants[46] = 0;
  constants[50] = 0;

  constants[51] = disp_num;
  wwritw(CHECKPOINTFILE,(char *) constants, MCONST*sizeof(int),100*sizeof(int),&junk);

  rclose(CHECKPOINTFILE,3);
  chkpt_close();
  free(constants);
  free(pointers);
  free(calcs);
  return;
}

