#include <stdio.h>
#include <stdlib.h>
#include <psio.h>
#include <libciomr/libciomr.h>
#include <file30.h>
#include <psifiles.h>

#include "input.h"
#include "defines.h"
#define EXTERN
#include "global.h"


void init_oldcalc()
{
  int s, ncalcs;
  
  /*--- initialize the old checkpoint file ---*/
  file30_init();

  /*--- check if MO vector is there ---*/
  ncalcs = file30_rd_ncalcs();
  if (ncalcs = 0) {
      fprintf(outfile,"  Did not find MOs in the checkpoint file\n");
      fprintf(outfile,"  Will continue without the MO projection\n");
      fprintf(stderr,"  Did not find MOs in the checkpoint file\n");
      fprintf(stderr,"  Will continue without the MO projection\n");
      chkpt_mos = 0;
      file30_close();
      return;
  }

  /*--- check if the symmetry is the same ---*/
  Oldcalc.symmetry = file30_rd_sym_label();
  if (strncmp(Oldcalc.symmetry,symmetry,3)) {
      fprintf(outfile,"  Old file30 uses different point group from the current one\n");
      fprintf(outfile,"  Will continue without the MO projection\n");
      fprintf(stderr,"  Old file30 uses different point group from the current one\n");
      fprintf(stderr,"  Will continue without the MO projection\n");
      chkpt_mos = 0;
      file30_close();
      return;
  }
  Oldcalc.nirreps = file30_rd_nirreps();

  /*--- read the old geometry ---*/
  Oldcalc.natom = file30_rd_natom();
  Oldcalc.geometry = file30_rd_geom();

  /*--- read the old basis ---*/
  Oldcalc.num_ao = file30_rd_nao();
  Oldcalc.num_so = file30_rd_nso();
  Oldcalc.sopi = file30_rd_sopi();
  Oldcalc.num_shells = file30_rd_nshell();
  Oldcalc.num_prims = file30_rd_nprim();
  Oldcalc.exponents = file30_rd_exps();
  Oldcalc.contr_coeff = file30_rd_contr_full();
  Oldcalc.shell_nucleus = file30_rd_snuc();
  Oldcalc.first_prim_shell = file30_rd_sprim();
  Oldcalc.shell_ang_mom = file30_rd_stype();
  Oldcalc.nprim_in_shell = file30_rd_snumg();
  Oldcalc.first_ao_shell = file30_rd_sloc();
  Oldcalc.usotao = file30_rd_usotao_new();

  /*--- Convert to C-numbering ---*/
  Oldcalc.max_angmom = 0;
  for(s=0;s<Oldcalc.num_shells;s++) {
      Oldcalc.shell_nucleus[s]--;
      Oldcalc.shell_ang_mom[s]--;
      Oldcalc.first_ao_shell[s]--;
      Oldcalc.first_prim_shell[s]--;
      if (Oldcalc.shell_ang_mom[s] > Oldcalc.max_angmom)
	  Oldcalc.max_angmom = Oldcalc.shell_ang_mom[s];
  }

  /*--- read the old SCF eigenvector, RHF case so far ---*/
  Oldcalc.ref = file30_rd_ref();
  if (Oldcalc.ref == ref_rhf   ||
      Oldcalc.ref == ref_rohf  ||
      Oldcalc.ref == ref_tcscf ||
      Oldcalc.ref == ref_rks)
      Oldcalc.spinrestr_ref = 1;
  else if (Oldcalc.ref == ref_uhf   ||
	   Oldcalc.ref == ref_uks)
      Oldcalc.spinrestr_ref = 0;
  else {
      fprintf(outfile,"  Input is not aware of reference type %d\n",Oldcalc.ref);
      fprintf(outfile,"  Will continue without the MO projection\n");
      chkpt_mos = 0;
      file30_close();
      return;
  }

  Oldcalc.iopen = file30_rd_iopen();
  Oldcalc.num_mo = file30_rd_nmo();
  Oldcalc.orbspi = file30_rd_orbspi();
  Oldcalc.clsdpi = file30_rd_clsdpi();
  Oldcalc.openpi = file30_rd_openpi();
  Oldcalc.escf = file30_rd_escf();
  if (Oldcalc.spinrestr_ref) { /* Spin-restricted case */
      Oldcalc.scf_evect_so = file30_rd_scf();
  }
  else { /* Spin-unrestricted case */
      Oldcalc.scf_evect_so_alpha = file30_rd_alpha_scf();
      Oldcalc.scf_evect_so_beta = file30_rd_beta_scf();
  }

  /*--- Done with the checkpoint file ---*/
  file30_close();
  
  return;
}


void cleanup_oldcalc()
{
  free_block(Oldcalc.geometry);
  free(Oldcalc.exponents);
  free_block(Oldcalc.contr_coeff);
  free(Oldcalc.shell_nucleus);
  free(Oldcalc.first_prim_shell);
  free(Oldcalc.shell_ang_mom);
  free(Oldcalc.nprim_in_shell);
  free(Oldcalc.first_ao_shell);
  free_block(Oldcalc.usotao);

  free(Oldcalc.orbspi);
  free(Oldcalc.clsdpi);
  free(Oldcalc.openpi);
  if (Oldcalc.spinrestr_ref)
      free(Oldcalc.scf_evect_so);
  else {
      free(Oldcalc.scf_evect_so_alpha);
      free(Oldcalc.scf_evect_so_beta);
  }

  return;
}


void save_oldmos()
{
  double **S12;

  if (max_angmom > Oldcalc.max_angmom)
    init_gto(max_angmom);
  else
    init_gto(Oldcalc.max_angmom);

  S12 = overlap_new_old();
  
  /*--- write things out ---*/
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_NEW);
  psio_write_entry(PSIF_OLD_CHKPT, "Old number of atoms", (char*) &Oldcalc.natom,
	     sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old geometry", (char*) Oldcalc.geometry[0],
	     Oldcalc.natom*3*sizeof(double));
  psio_write_entry(PSIF_OLD_CHKPT, "Old number of irreps", (char*) &Oldcalc.nirreps,
	     sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old sopi", (char*) Oldcalc.sopi,
	     Oldcalc.nirreps*sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old reference type", (char*) &Oldcalc.ref,
	     sizeof(reftype));
  psio_write_entry(PSIF_OLD_CHKPT, "Old num_so", (char*) &Oldcalc.num_so,
	     sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old num_mo", (char*) &Oldcalc.num_mo,
	     sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old orbspi", (char*) Oldcalc.orbspi,
	     Oldcalc.nirreps*sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old clsdpi", (char*) Oldcalc.clsdpi,
	     Oldcalc.nirreps*sizeof(int));
  psio_write_entry(PSIF_OLD_CHKPT, "Old openpi", (char*) Oldcalc.openpi,
	     Oldcalc.nirreps*sizeof(int));
  if (Oldcalc.spinrestr_ref) { /* Spin-restricted case */
    psio_write_entry(PSIF_OLD_CHKPT, "Old MOs", (char*) Oldcalc.scf_evect_so[0],
	       Oldcalc.num_so*Oldcalc.num_mo*sizeof(double));
  }
  else { /* Spin-unrestricted case */
    psio_write_entry(PSIF_OLD_CHKPT, "Old alpha MOs", (char*) Oldcalc.scf_evect_so_alpha[0],
	       Oldcalc.num_so*Oldcalc.num_mo*sizeof(double));
    psio_write_entry(PSIF_OLD_CHKPT, "Old beta MOs", (char*) Oldcalc.scf_evect_so_alpha[0],
	       Oldcalc.num_so*Oldcalc.num_mo*sizeof(double));
  }
  psio_write_entry(PSIF_OLD_CHKPT, "New-Old overlap", (char*) S12[0],
	     num_so*Oldcalc.num_so*sizeof(double));
  psio_close(PSIF_OLD_CHKPT, 1);

  free_block(S12);
}
