#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <physconst.h>
#include <psifiles.h>
#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"
#include "globals.h"

void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_params(void);
void init_ioff(void);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
struct dpd_file4_cache_entry *priority_list(void);
void energy(void);
void opdm(void);
void cleanup(void);
void exit_io(void);

int main(int argc, char *argv[])
{
  int *cachefiles;
  int **cachelist;
 
  struct dpd_file4_cache_entry *priority;
  
  init_io(argc,argv);
  title();

  get_moinfo();
  get_params();
  init_ioff();
  
  cachefiles = init_int_array(PSIO_MAXUNIT);

  cachelist = cacheprep_rhf(params.cachelev, cachefiles);

  priority = priority_list();
      
  dpd_init(0, mo.nirreps, params.memory, params.cachetype, cachefiles,
           cachelist, priority, 2, mo.actdoccpi, mo.actdoccsym,
	   mo.actvirtpi, mo.actvirtsym);
  
  energy();
  
  if (params.opdm) {
  opdm();
  }
  
  dpd_close(0);
  
  cleanup();

  exit_io();
  
  exit(0);
}

void init_io(int argc, char *argv[])
{
  int i;
  extern char *gprgid();
  char *progid;

  progid = (char *)malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  
  psio_init();
  for(i=CC_MIN; i <= CC_LAMPS; i++) 
    psio_open(i,1);
  for(i=CC_HBAR; i <= CC_MAX; i++) 
    psio_open(i,0);
}

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*          MP2          *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
}

/* Variable Representation
 * n  -> number of 
 * pi -> per irrep
 * fz -> frozen
 * act -> active
 * docc -> doubly occupied
 * virt -> virtual
 * mo -> molecular orbitals
 * when all docc or all virt are used the
 * meaning is intrinsic
*/

void get_moinfo(void)
{
  int i;
  
  chkpt_init(PSIO_OPEN_OLD);

  mo.nmo = chkpt_rd_nmo();
  mo.nso = chkpt_rd_nso();
  mo.nirreps = chkpt_rd_nirreps();
  mo.irreplabels = chkpt_rd_irr_labs();

  mo.mopi = chkpt_rd_orbspi();
  mo.doccpi = chkpt_rd_clsdpi();
  
  mo.Enuc = chkpt_rd_enuc();
  mo.Escf = chkpt_rd_escf();
  mo.scfevals = chkpt_rd_evals();
  
  mo.fzdoccpi = get_frzcpi();
  mo.fzvirtpi = get_frzvpi();
  
  chkpt_close();

  mo.nfzdocc = 0;
  mo.nfzvirt = 0;
  for (i=0; i<mo.nirreps; i++) {
    mo.nfzdocc += mo.fzdoccpi[i];
    mo.nfzvirt += mo.fzvirtpi[i];
  }

  mo.virtpi = init_int_array(mo.nirreps);
  for(i=0; i < mo.nirreps; i++) {
    mo.virtpi[i] = mo.mopi[i]-mo.doccpi[i];
  }
  
  mo.ndocc = 0;
  mo.nvirt = 0;
  for(i=0; i < mo.nirreps; i++) {
    mo.ndocc += mo.doccpi[i];
    mo.nvirt += mo.mopi[i] - mo.doccpi[i];
  }
  
  mo.actdoccpi = init_int_array(mo.nirreps);
  mo.actvirtpi = init_int_array(mo.nirreps);
  for(i=0; i < mo.nirreps; i++) {
    mo.actdoccpi[i] = mo.doccpi[i]-mo.fzdoccpi[i];
    mo.actvirtpi[i] = mo.virtpi[i]-mo.fzvirtpi[i];
  }

  mo.nactdocc = 0;
  mo.nactvirt = 0;
  for (i=0; i < mo.nirreps; i++) {
    mo.nactdocc += mo.actdoccpi[i];
    mo.nactvirt += mo.actvirtpi[i];
  }
 
  mo.nactmo = mo.nactdocc + mo.nactvirt;
  
  mo.actdoccsym = init_int_array(mo.nactmo);
  mo.actvirtsym = init_int_array(mo.nactmo);
  psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
	         (char *) mo.actdoccsym, sizeof(int)*mo.nactmo);
  psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
		 (char *) mo.actvirtsym, sizeof(int)*mo.nactmo);

  mo.docc_off = init_int_array(mo.nirreps);
  mo.virt_off = init_int_array(mo.nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orb Offsets",
	          (char *) mo.docc_off, sizeof(int)*mo.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orb Offsets",
		  (char *) mo.virt_off, sizeof(int)*mo.nirreps);

  mo.qt_docc = init_int_array(mo.nactmo);
  mo.qt_virt = init_int_array(mo.nactmo);

  psio_read_entry(CC_INFO, "CC->QT Active Occ Order",
                 (char *) mo.qt_docc, sizeof(int)*mo.nactmo);
  psio_read_entry(CC_INFO, "CC->QT Active Virt Order",
	         (char *) mo.qt_virt, sizeof(int)*mo.nactmo);
	      
  fprintf(outfile,"\n");
  fprintf(outfile,"\tChkpt Parameters:\n");
  fprintf(outfile,"\t--------------------\n");
  fprintf(outfile,"\tNumber of irreps     = %d\n",mo.nirreps);
  fprintf(outfile,"\tNumber of MOs        = %d\n",mo.nmo);
  fprintf(outfile,"\n");
  fprintf(outfile,
    "\tLabel\tFZDC\tACTD\tDOCC\tACTV\tFZVI\tVIRT\tMOs\n");
  fprintf(outfile,
    "\t-----\t----\t----\t----\t----\t----\t----\t---\n");
  for(i=0; i < mo.nirreps; i++) {
    fprintf(outfile,
    "\t  %s \t  %d\t  %d\t  %d\t  %d\t  %d\t  %d\t %d\n",
	    mo.irreplabels[i],mo.fzdoccpi[i],mo.actdoccpi[i],mo.doccpi[i],
	    mo.actvirtpi[i],mo.fzvirtpi[i],mo.virtpi[i],mo.mopi[i]);
  }
  
  fprintf(outfile,"\n");
  fprintf(outfile,"\tNuclear Rep. energy \t=\t  %.12f\n",mo.Enuc);
  fprintf(outfile,"\tSCF energy          \t=\t%.12f\n",mo.Escf);
}


void get_params()
{
  int errcod;
  char *cachetype = NULL;
  
  errcod = ip_string("WFN", &(params.wfn), 0);

  errcod = ip_string("REFERENCE", &(params.ref),0);
  if (strcmp(params.ref,"RHF")) {
    fprintf(outfile, "\nIncorrect Reference: RHF only\n");
    abort();
  }
  
  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);
  
  params.opdm = 0;
  errcod = ip_boolean("OPDM", &(params.opdm),0);
  
  if (params.opdm) {
    params.opdm_write = 1;
  }  
  else {
    params.opdm_write = 0;
  }
  errcod = ip_boolean("OPDM_WRITE", &(params.opdm_write),0);
  
  params.opdm_print = 0;
  errcod = ip_boolean("OPDM_PRINT", &(params.opdm_print),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  
  params.cachetype = 1;
  errcod = ip_string("CACHETYPE", &(cachetype),0);
  if (cachetype != NULL && strlen(cachetype)) {
    if (!strcmp(cachetype,"LOW")) 
      params.cachetype = 1;
    else if (!strcmp(cachetype,"LRU")) 
      params.cachetype = 0;
    else {
      fprintf(outfile, "Invalide CACHETYPE = %s\n",cachetype);
      abort();
    }
    free(cachetype);
  }
  
  fndcor(&(params.memory),infile,outfile);
 
  fprintf(outfile, "\n");
  fprintf(outfile, "\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function \t=\t%s\n", params.wfn);
  fprintf(outfile, "\tReference WFN \t=\t%s\n", params.ref);
  fprintf(outfile, "\tCache Level   \t=\t%d\n", params.cachelev);
  fprintf(outfile, "\tCache Type    \t=\t%s\n", params.cachetype ? "LOW":"LRU");
  fprintf(outfile, "\tMemory (MB)   \t=\t%.1f\n",params.memory/1e6);
  fprintf(outfile, "\tPrint Level   \t=\t%d\n", params.print);
  fprintf(outfile, "\tOPDM          \t=\t%s\n", params.opdm ? "YES":"NO");
  fprintf(outfile, "\tWrite OPDM    \t=\t%s\n", params.opdm_write ? "YES":"NO");
  fprintf(outfile, "\tPrint OPDM    \t=\t%s\n", params.opdm_print ? "YES":"NO");
}

void init_ioff(void)
{
  int i;
  
  mo.ioff = init_int_array(MAXIOFF);
  mo.ioff[0] = 0;
  for(i=1; i < MAXIOFF; i++) {
    mo.ioff[i] = mo.ioff[i-1] + i;
  }
  
}

void cleanup(void)
{
  int i;
  
  free(params.wfn);
  free(params.ref);
  free(mo.fzdoccpi);
  free(mo.fzvirtpi);
  free(mo.actdoccpi);
  free(mo.actdoccsym);
  free(mo.actvirtpi);
  free(mo.actvirtsym);
  free(mo.doccpi);
  free(mo.virtpi);
  free(mo.mopi);
  for(i=0; i < mo.nirreps; i++)
    free(mo.irreplabels[i]);
  free(mo.irreplabels);
  free(mo.ioff);
  free(mo.scfevals);
}

void exit_io(void)
{
  int i;

  for(i=CC_MIN; i <= CC_MAX; i++) 
    psio_close(i,1);
  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
  char *prgid = "MP2";
  return(prgid);
}
