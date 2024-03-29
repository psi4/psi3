/*! \defgroup RESPONSE response: Compute various response properties */

/*! 
** \file
** \ingroup RESPONSE
** \brief Module to compute various response properties
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

namespace psi { namespace response {

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io(int argc, char *argv[]);
void init_ioff(void);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void exit_io(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void build_A_RHF(void);
void build_B_RHF(void);
void polar(void);
void optrot(void);
void roa(void);
void transpert(const char *pert);

}} // namespace psi::response

using namespace psi::response;

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;

  init_io(argc, argv);
  init_ioff();
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, 
	     NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi, moinfo.avir_sym,
	     moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  } 
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  }

  if(params.ref == 0) {

    build_A_RHF();
    build_B_RHF();
    
    if(!strcmp(params.prop,"POLARIZABILITY")) polar();
    else if(!strcmp(params.prop,"ROTATION")) optrot();
    else if(!strcmp(params.prop,"ROA")) roa();
  }

  dpd_close(0);
  cleanup();
  exit_io();
  exit(0);
}

extern "C" {const char *gprgid() { const char *prgid = "RESPONSE"; return(prgid); }}

namespace psi { namespace response {

void init_io(int argc, char *argv[])
{
  int i;
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0); /* this assumes no cmdline args except filenames */
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  psio_init(); psio_ipv1_config();

  psio_open(PSIF_MO_HESS,0);
  psio_open(CC_INFO, PSIO_OPEN_OLD);
  psio_open(CC_OEI, PSIO_OPEN_OLD);
  psio_open(CC_CINTS, PSIO_OPEN_OLD);
  psio_open(CC_DINTS, PSIO_OPEN_OLD);
  psio_open(CC_TMP0, PSIO_OPEN_NEW);
}

void exit_io(void)
{
  int i;
 
  psio_close(PSIF_MO_HESS,1);
  psio_close(CC_INFO, 1);
  psio_close(CC_OEI, 1);
  psio_close(CC_CINTS, 1);
  psio_close(CC_DINTS, 1);
  psio_close(CC_TMP0, 1);

  psio_done();
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

}} // namespace psi::response
