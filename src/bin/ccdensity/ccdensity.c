/*
**  CCDENSITY: Program to calculate the coupled-cluster one- and 
**             two-particle densities.
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <math.h>
#include <psifiles.h>
#include "globals.h"

/* Function prototypes */
void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_frozen(void);
void get_params(void);
void exit_io(void);
void onepdm(void);
void sortone(void);
void twopdm(void);
void energy(void);
void resort_tei(void);
void resort_gamma(void);
void lag(void);
void build_X(void);
void build_A(void);
void build_Z(void);
void relax_I(void);
void relax_D(void);
void sortI(void);
void fold(void);
void deanti(void);
void add_ref_ROHF(struct iwlbuf *);
void add_ref_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void add_core_ROHF(struct iwlbuf *);
void add_core_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void dump_ROHF(struct iwlbuf *);
void dump_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void kinetic(void);
void dipole(void);
void probable(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void add_eom_d(void);
void sort_Ls(void);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  struct iwlbuf OutBuf;
  struct iwlbuf OutBuf_AA, OutBuf_BB, OutBuf_AB;
  
  init_io(argc,argv);
  title();
  get_moinfo();
  /*  get_frozen(); */
  get_params();

  if(moinfo.nfzc || moinfo.nfzv) {
    fprintf(outfile, "\n\tGradients involving frozen orbitals not yet available.\n");
    exit(PSI_RETURN_FAILURE);
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

  }
  else if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
	     cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
	     moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  if (params.ground)
    sort_Ls();

  onepdm();

  if (!params.ground) {
    add_eom_d(); /* R0 * ground D + excited parts of onepdm */
    fprintf(outfile,"\tStill using ground-state 2-pdm\n");
  }


  twopdm();

  if(!params.aobasis) energy();
  sortone();
  /* dipole(); */
  kinetic();

  /*
  fprintf(outfile,"\tDipole moments without orbital relaxation\n");
  dipole();
  */

  lag();

  /*

  dpd_init(1, moinfo.nirreps, params.memory, 2, frozen.occpi, frozen.occ_sym,
  frozen.virtpi, frozen.vir_sym);
  */

  /*  if(moinfo.nfzc || moinfo.nfzv) {
      resort_gamma();
      resort_tei();
      } */

  build_X();
  build_A();
  build_Z();

  relax_I();
  if(params.relax_opdm) {
    relax_D();
  }
  sortone();

  sortI();

  fold();
  deanti();

  /*  
      dpd_close(0);
      dpd_close(1);
  */

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

    add_core_ROHF(&OutBuf);
    add_ref_ROHF(&OutBuf);
    dump_ROHF(&OutBuf);

    iwl_buf_flush(&OutBuf, 1);
    iwl_buf_close(&OutBuf, 1);

  }
  else if(params.ref == 2) { /** UHF **/

    iwl_buf_init(&OutBuf_AA, PSIF_MO_AA_TPDM, params.tolerance, 0, 0);
    iwl_buf_init(&OutBuf_BB, PSIF_MO_BB_TPDM, params.tolerance, 0, 0);
    iwl_buf_init(&OutBuf_AB, PSIF_MO_AB_TPDM, params.tolerance, 0, 0);

    /*    add_core_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB); */
    add_ref_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB);
    dump_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB);

    iwl_buf_flush(&OutBuf_AA, 1);
    iwl_buf_flush(&OutBuf_BB, 1);
    iwl_buf_flush(&OutBuf_AB, 1);
    iwl_buf_close(&OutBuf_AA, 1);
    iwl_buf_close(&OutBuf_BB, 1);
    iwl_buf_close(&OutBuf_AB, 1);

  }

  if(params.ref == 2) { 
    dpd_close(0);
    cleanup();
    exit_io();
    exit(PSI_RETURN_SUCCESS);
  }

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup(); 
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

void init_io(int argc, char *argv[])
{
  int i, num_unparsed;
  extern char *gprgid();
  char *progid, *argv_unparsed[100];;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  params.ground = 1;
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--excited")) {
      params.ground = 0;
    }
    else {
      argv_unparsed[num_unparsed++] = argv[i];
    }
  }

  psi_start(num_unparsed,argv_unparsed,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  psio_init();

  /* Open all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);

  if (!params.ground) {
    /* assume symmetry of L is that of R */
    psio_read_entry(CC_INFO,"EOM R0", (char *) &(params.R0),sizeof(double));
    psio_read_entry(CC_INFO,"CCEOM Energy", (char *) &(params.cceom_energy),sizeof(double));
  }
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCDENSITY       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;
 
  /* Close all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);

  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "CCDENSITY";

   return(prgid);
}

/* R0 * ground D + excited parts of onepdm */
void add_eom_d(void) {
  dpdfile2 DG, DX;

  fprintf(outfile,"Multiplying density by R0 and adding in excited parts of onepdm\n");

  dpd_file2_init(&DG, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 0, 0, "DIJ");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  dpd_file2_init(&DG, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 0, 0, "Dij");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  dpd_file2_init(&DG, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 1, 1, "DAB");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  dpd_file2_init(&DG, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 1, 1, "Dab");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  dpd_file2_init(&DG, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 0, 1, "DIA");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  dpd_file2_init(&DG, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 0, 1, "Dia");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  dpd_file2_init(&DG, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 0, 1, "DAI");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);
  
  dpd_file2_init(&DG, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_scm(&DG, params.R0);
  dpd_file2_init(&DX, EOM_D, 0, 0, 1, "Dai");
  dpd_file2_axpy(&DX, &DG, 1.0, 0);
  dpd_file2_close(&DX);
  dpd_file2_close(&DG);

  return;
}
