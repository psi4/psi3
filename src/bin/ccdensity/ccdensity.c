/*
**  CCDENSITY: Program to calculate the coupled-cluster one- and 
**             two-particle densities.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
void setup_LR(void);
void G_build(void);
void x_oe_intermediates(void);
void x_onepdm(void);
void x_te_intermediates(void);
void x_Gijkl(void);
void x_Gabcd(void);
void x_Gibja(void);
void x_Gijka(void);
void x_Gijab_ROHF(void);
void x_Gciab(void);
void V_build_x(void);
void x_xi1(void);
void x_xi_zero(void);
void x_xi2(void);
void x_xi_oe_intermediates(void);
void G_norm(void);
void zero_onepdm(void);
void zero_twopdm(void);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  struct iwlbuf OutBuf;
  struct iwlbuf OutBuf_AA, OutBuf_BB, OutBuf_AB;
  dpdfile2 D;
  
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

  /* CC_GLG will contain L, or R0*L + Zeta, if relaxed and zeta is available */
  /* CC_GL will contain L */
  setup_LR();

  /* compute ground state parts of onepdm or put zeroes there */
  if ( (!params.calc_xi) && ( (params.L_irr == params.G_irr) || (params.use_zeta) ) )
    onepdm();
  else
    zero_onepdm();

  if ( (!params.calc_xi) && ( (params.L_irr == params.G_irr) || (params.use_zeta) ) ) {
    V_build(); /* uses CC_GLG, writes to CC_MISC */
    G_build(); /* uses CC_GLG, writes to CC_GLG */
  }

  /* calculate intermediates not already on disk */
  if (!params.restart) {
    if (!params.ground) {
      /* the following intermediates go in EOM_TMP */
      V_build_x(); /* use CC_GL write to EOM_TMP */
      x_oe_intermediates();
      x_te_intermediates();
      if (params.calc_xi) {
        /* the following intermediates go in EOM_TMP_XI */
        x_xi_intermediates();
      }
    }
  }

  /* If calculating xi, do it and then quit */
  /* leave intermediates in CC_MISC and EOM_TMP on disk */
  if ( params.calc_xi ) {
    x_xi_zero(); /* make blank Xi to help debugging, in x_xi1.c */
    x_xi1();
    x_xi2();
    dpd_close(0);
    if(params.ref == 2) cachedone_uhf(cachelist);
    else cachedone_rhf(cachelist);
    free(cachefiles);
    cleanup();
    psio_close(EOM_TMP_XI,0);
    psio_open(EOM_TMP_XI,PSIO_OPEN_NEW);
    exit_io();
    exit(PSI_RETURN_SUCCESS);
  }

  if ( (!params.calc_xi) && ( (params.L_irr == params.G_irr) || (params.use_zeta) ) )
    twopdm();
  else
    zero_twopdm();

  //fprintf(outfile,"After ground state parts\n");
  //G_norm();

  /* add in non-R0 parts of onepdm and twopdm */
  if (!params.ground) {
    x_onepdm();
    x_Gijkl();
    x_Gabcd();
    x_Gibja();
    x_Gijka();
    x_Gciab();
    x_Gijab_ROHF();
  }
  //fprintf(outfile,"After excited state parts\n");
  //G_norm();
  
  if(!params.aobasis) energy();

  sortone();
  // dipole();
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
  /*testing
  fprintf(outfile,"After orbital response\n");
  if(!params.aobasis) energy();
  */

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
  params.calc_xi = 0;
  params.restart = 0;
  params.use_zeta = 0;
  params.user_transition = 0;
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--excited")) {
      params.ground = 0;
    }
    else if (!strcmp(argv[i],"--use_zeta")) {
      params.use_zeta = 1;
      params.ground = 0;
      params.restart = 1;
    }
    else if (!strcmp(argv[i],"--calc_xi")) {
      params.calc_xi = 1;
      params.ground = 0;
      params.restart = 0;
    }
    else if (!strcmp(argv[i],"--transition")) {
      params.user_transition = 1;
      sscanf(argv[++i], "%d",&(params.L_irr));
      sscanf(argv[++i], "%d",&(params.L_root));
      sscanf(argv[++i], "%d",&(params.R_root));
      params.R_irr = params.L_irr; /* assume same for now */
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
  /* erase files for easy debugging */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,PSIO_OPEN_OLD);
  /*
  psio_close(CC_GR,0);
  psio_close(CC_GL,0);
  psio_close(EOM_TMP0,0);
  psio_open(CC_GR,PSIO_OPEN_NEW);
  psio_open(CC_GL,PSIO_OPEN_NEW);
  psio_open(EOM_TMP0,PSIO_OPEN_NEW);
  */
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

  /* delete temporary EOM files */
  psio_close(EOM_TMP0,0);
  psio_close(EOM_TMP1,0);
//psio_close(CC_GLG,0);
  psio_open(EOM_TMP0,PSIO_OPEN_NEW);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
//psio_open(CC_GLG,PSIO_OPEN_NEW);
  if (!params.calc_xi) {
    psio_close(EOM_TMP,0);
    psio_open(EOM_TMP,PSIO_OPEN_NEW);
  }

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

