/*
** TRANSQT:
** Program to transform one- and two-electron integrals over
** symmetry-orbitals to integrals over molecular orbitals.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "globals.h"

/* Function prototypes */
void init_io(int argc, char *argv[]);
void init_ioff(void);
void title(void);
void get_params(void);
void get_moinfo(void);
void cleanup(void);
void exit_io(void);

int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhd(int **cachelist);

int file_build_presort(dpdfile4 *, int, double, long int, int, 
		       int, double **, double **, double **, double **, int);

void transform_two_rhf(void);
#define PSIF_HALFT0 91
#define PSIF_HALFT1 92

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  dpdfile4 I;
  dpdbuf4 J, K;
  int h, pq, rs, p, q, r, s, Gs, Gr, PQ, RS, i;
  int nrows, ncols, nlinks;
  double **TMP, **H, **D, **F, *oei;
  double ***C;
  struct iwlbuf MBuff;
  int stat, noei;
  int *offset;
  double efzc;
  double **scratch;

  init_io(argc,argv);
  init_ioff();
  title();
  get_params();
  get_moinfo();

  TMP = block_matrix(moinfo.nso,moinfo.nso);

  cachefiles = init_int_array(PSIO_MAXUNIT);
  if(params.ref == 2) cachelist = cacheprep_uhf(params.cachelev, cachefiles);
  else cachelist = cacheprep_rhf(params.cachelev, cachefiles);

  dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist,
	   NULL, 2, moinfo.sopi, moinfo.sosym, moinfo.mopi, moinfo.mosym);

  psio_open(PSIF_SO_PRESORT, 0);
  psio_open(PSIF_HALFT0, 0);
  psio_open(PSIF_HALFT1, 0);

  /*** Starting one-electron transforms ***/

  /* For the one-electron integral transform, the full MO list is best */
  C = (double ***) malloc(1 * sizeof(double **));
  chkpt_init(PSIO_OPEN_OLD);
  C[0] = chkpt_rd_scf();
  chkpt_close();

  /* build the frozen-core density (RHF) */
  offset = init_int_array(moinfo.nirreps);
  for(h=1; h < moinfo.nirreps; h++)
    offset[h] = offset[h-1] + moinfo.sopi[h-1];

  D = block_matrix(moinfo.nso,moinfo.nso);
  for(h=0; h < moinfo.nirreps; h++)
    for(p=offset[h]; p < offset[h]+moinfo.sopi[h]; p++)
      for(q=offset[h]; q < offset[h]+moinfo.sopi[h]; q++)
	for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++)
	  D[p][q] += C[0][p][i] * C[0][q][i];

  free(offset);

  /* pre-sort the SO-basis two-electron integrals and generate the fzc operator */
  F = block_matrix(moinfo.nso,moinfo.nso);
  fprintf(outfile, "\n\tPresorting SO-basis two-electron integrals.\n");
  fflush(outfile);
  dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (pq,rs)");
  file_build_presort(&I, PSIF_SO_TEI, params.tolerance, params.memory, 1, 
		     moinfo.nfzc, D, NULL, F, NULL, 0);
  dpd_file4_close(&I);

  /* read the bare one-electron integrals */
  noei = moinfo.nso*(moinfo.nso+1)/2;
  oei = init_array(noei);
  H = block_matrix(moinfo.nso,moinfo.nso);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_T, oei, noei, 0, 0, outfile);
  for(p=0,pq=0; p < moinfo.nso; p++)
    for(q=0; q <= p; q++,pq++)
      H[p][q] = H[q][p] = oei[pq];
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_V, oei, noei, 0, 0, outfile);
  for(p=0,pq=0; p < moinfo.nso; p++)
    for(q=0; q <= p; q++,pq++) {
      H[p][q] += oei[pq];
      if(p!=q) H[q][p] += oei[pq];
    }

  /* add the remaining one-electron terms to the fzc operator */
  for(p=0; p < moinfo.nso; p++)
    for(q=0; q < moinfo.nso; q++)
      F[p][q] += H[p][q];

  /* compute the frozen-core energy and write it to the chkpt file*/
  efzc = 0.0;
  if(params.ref == 0 || params.ref == 1) {
    for(p=0; p < moinfo.nso; p++)
      for(q=0; q < moinfo.nso; q++)
	efzc += D[p][q] * (H[p][q] + F[p][q]);
  }
  fprintf(outfile, "\tFrozen-core energy = %20.15f\n", efzc);
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_efzc(efzc);
  chkpt_close();

  /* transform the bare one-electron integrals */
  nrows = moinfo.nso;
  ncols = moinfo.nmo;
  nlinks = moinfo.nso;
  if(nrows && ncols && nlinks)
    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,H[0],nlinks,C[0][0],ncols,
	    0.0,TMP[0],moinfo.nso);

  nrows = moinfo.nmo;
  ncols = moinfo.nmo;
  nlinks = moinfo.nso;
  if(nrows && ncols && nlinks)
    C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[0][0],ncols,TMP[0],moinfo.nso,
	    0.0,H[0],nlinks);

  if(params.print_lvl > 2) {
    fprintf(outfile, "\tOne-electron integrals (MO basis):\n");
    mat_print(H, moinfo.nmo, moinfo.nmo, outfile);
  }

  for(p=0; p < moinfo.nmo; p++)
    for(q=0; q < moinfo.nmo; q++) {
      pq = INDEX(moinfo.pitzer2qt[p],moinfo.pitzer2qt[q]);
      oei[pq] = H[p][q];
    }

  noei = moinfo.nmo*(moinfo.nmo+1)/2;
  iwl_wrtone(PSIF_OEI, PSIF_MO_OEI, noei, oei);

  /* transform the frozen-core operator */
  nrows = moinfo.nso;
  ncols = moinfo.nmo;
  nlinks = moinfo.nso;
  if(nrows && ncols && nlinks)
    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,F[0],nlinks,C[0][0],ncols,
	    0.0,TMP[0],moinfo.nso);

  nrows = moinfo.nmo;
  ncols = moinfo.nmo;
  nlinks = moinfo.nso;
  if(nrows && ncols && nlinks)
    C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[0][0],ncols,TMP[0],moinfo.nso,
	    0.0,F[0],nlinks);

  if(params.print_lvl > 1) {
    fprintf(outfile, "\tFrozen-core operator (MO basis):\n");
    mat_print(F, moinfo.nmo, moinfo.nmo, outfile);
  }

  for(p=0; p < moinfo.nmo; p++)
    for(q=0; q < moinfo.nmo; q++) {
      pq = INDEX(moinfo.pitzer2qt[p],moinfo.pitzer2qt[q]);
      oei[pq] = F[p][q];
    }

  iwl_wrtone(PSIF_OEI, PSIF_MO_FZC, noei, oei);

  free(oei);
  free_block(H);
  free_block(F);
  free_block(D);

  free_block(C[0]);
  free(C);

  /*** One-electron transforms complete ***/

  /*** Starting two-electron transforms ***/

  /* For the two-electron integral transform, symmetry-blocked MOs are best */
  C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  chkpt_init(PSIO_OPEN_OLD);
  for(h=0; h < moinfo.nirreps; h++) {
    scratch = chkpt_rd_scf_irrep(h);
    C[h] = block_matrix(moinfo.sopi[h],moinfo.mopi[h]);
    for(q=0; q < moinfo.mopi[h]; q++)
      for(p=0; p < moinfo.sopi[h]; p++)
	C[h][p][q] = scratch[p][q+moinfo.frdocc[h]];
    if(params.print_lvl > 2) {
      fprintf(outfile, "\tMOs for irrep %d:\n",h);
      mat_print(C[h], moinfo.sopi[h], moinfo.mopi[h], outfile);
    }
    free_block(scratch);
  }
  chkpt_close();

  fprintf(outfile, "\tStarting first half-transformation.\n");
  fflush(outfile);
  dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_row_init(&J, h);
    dpd_buf4_mat_irrep_row_init(&K, h);

    for(pq=0; pq < J.params->rowtot[h]; pq++) {
      dpd_buf4_mat_irrep_row_rd(&J, h, pq);

      for(Gr=0; Gr < moinfo.nirreps; Gr++) {
	Gs = h^Gr;

	nrows = moinfo.sopi[Gr];
	ncols = moinfo.mopi[Gs];
	nlinks = moinfo.sopi[Gs];
	rs = J.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][0][rs],nlinks,
		  C[Gs][0],ncols,0.0,TMP[0],moinfo.nso);

	nrows = moinfo.mopi[Gr];
	ncols = moinfo.mopi[Gs];
	nlinks = moinfo.sopi[Gr];
	rs = K.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nrows,TMP[0],moinfo.nso,
		  0.0,&K.matrix[h][0][rs],ncols);
      }

      dpd_buf4_mat_irrep_row_wrt(&K, h, pq);
    }

    dpd_buf4_mat_irrep_row_close(&J, h);
    dpd_buf4_mat_irrep_row_close(&K, h);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  psio_close(PSIF_SO_PRESORT, 0);

  fprintf(outfile, "\tSorting half-transformed integrals.\n");
  fflush(outfile);
  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 8, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, PSIF_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  psio_close(PSIF_HALFT0, 0);

  fprintf(outfile, "\tStarting second half-transformation.\n");
  fflush(outfile);
  iwl_buf_init(&MBuff, PSIF_MO_TEI, params.tolerance, 0, 0);

  dpd_buf4_init(&J, PSIF_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_row_init(&J, h);
    dpd_buf4_mat_irrep_row_init(&K, h);

    for(pq=0; pq < J.params->rowtot[h]; pq++) {
      dpd_buf4_mat_irrep_row_rd(&J, h, pq);

      for(Gr=0; Gr < moinfo.nirreps; Gr++) {
	Gs = h^Gr;

	nrows = moinfo.sopi[Gr];
	ncols = moinfo.mopi[Gs];
	nlinks = moinfo.sopi[Gs];
	rs = J.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][0][rs],nlinks,
		  C[Gs][0],ncols,0.0,TMP[0],moinfo.nso);

	nrows = moinfo.mopi[Gr];
	ncols = moinfo.mopi[Gs];
	nlinks = moinfo.sopi[Gr];
	rs = K.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nrows,TMP[0],moinfo.nso,
		  0.0,&K.matrix[h][0][rs],ncols);
      }

      p = moinfo.act2qt[K.params->roworb[h][pq][0]];
      q = moinfo.act2qt[K.params->roworb[h][pq][1]];
      PQ = INDEX(p,q);
      for(rs=0; rs < K.params->coltot[h]; rs++) {
	r = moinfo.act2qt[K.params->colorb[h][rs][0]];
	s = moinfo.act2qt[K.params->colorb[h][rs][1]];
	RS = INDEX(r,s);
	if(r >= s && RS <= PQ)
	  iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][0][rs], (params.print_lvl>10), outfile, 0);
      }

      /*       dpd_buf4_mat_irrep_row_wrt(&K, h, pq); */
    }

    dpd_buf4_mat_irrep_row_close(&J, h);
    dpd_buf4_mat_irrep_row_close(&K, h);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  for(h=0; h < moinfo.nirreps; h++)
    free_block(C[h]);
  free(C);

  psio_close(PSIF_HALFT1, 0);

  fprintf(outfile, "\tTwo-electron integral transformation complete.\n");
  fflush(outfile);

  /*** Two-electron transforms complete ***/

  dpd_close(0);

  free_block(TMP);

  cleanup();
  free(ioff);
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

char *gprgid(void) { char *prgid = "TRANSQT"; return (prgid); }

void init_io(int argc, char *argv[])
{
  int i;
  extern char *gprgid(void);
  char *progid;

  psi_start(argc-1,argv+1,0);

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s", gprgid());
  ip_cwk_add(progid);
  free(progid);

  tstart(outfile);
  psio_init();

  psio_open(CC_INFO, PSIO_OPEN_NEW);
  psio_open(CC_MISC, PSIO_OPEN_NEW);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile,"\t**************************************************\n");
  fprintf(outfile,"\t* TRANSQT:  Program to transform integrals from  *\n");
  fprintf(outfile,"\t*           the SO basis to the MO basis.        *\n");
  fprintf(outfile,"\t*                                                *\n");
  fprintf(outfile,"\t*            Daniel, David, & Justin             *\n");
  fprintf(outfile,"\t**************************************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  psio_close(CC_INFO,1);
  psio_close(CC_MISC,1);
  psio_done();
  tstop(outfile);
  psi_stop();
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  for(i=1; i < IOFF_MAX; i++)
    ioff[i] = ioff[i-1] + i;
}
