#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>

FILE *infile, *outfile;
char *psi_file_prefix;

void init_io(int argc, char *argv[]);
void exit_io(void);

int main(int argc, char *argv[])
{
  int i, j, h, I;
  int this_atom;
  FILE *molfile;
  int natom, nprim, nshell;
  int nmo, nso, nao, nirreps, ndocc;
  int *snumg, *stype, *snuc, *sprim;
  int *orbspi, *clsdpi, *orbsym;
  double **scf_so, **scf;
  double *evals;
  double **usotao;
  double **geom, *zvals, *exps, *contr, occ;
  char **atom_labs, *am_labs[7], **irr_labs;

  /* I don't think molden knows anything > f, but it's here anyway */
  am_labs[0] = strdup("s");
  am_labs[1] = strdup("p");
  am_labs[2] = strdup("d");
  am_labs[3] = strdup("f");
  am_labs[4] = strdup("g");
  am_labs[5] = strdup("h");
  am_labs[6] = strdup("i");

  init_io(argc, argv);

  ffile(&molfile, "molden.in", 0);
  fprintf(molfile, "[Molden Format]\n");

  chkpt_init(PSIO_OPEN_OLD);

  /** coordinates **/
  fprintf(molfile, "[Atoms] AU\n");

  natom = chkpt_rd_natom();
  geom = chkpt_rd_geom();
  zvals = chkpt_rd_zvals();
  atom_labs = chkpt_rd_felement();

  for(i=0; i < natom; i++)
    fprintf(molfile, "%s %3d %3d %20.10f %20.10f %20.10f\n", atom_labs[i], i+1, (int) zvals[i],
	    geom[i][0], geom[i][1], geom[i][2]);

  free_block(geom);
  free(zvals);
  for(i=0; i < natom; i++) free(atom_labs[i]);
  free(atom_labs);

  /** basis set **/
  fprintf(molfile, "[GTO]\n");

  nprim = chkpt_rd_nprim();
  nshell = chkpt_rd_nshell();
  snumg = chkpt_rd_snumg();
  exps = chkpt_rd_exps();
  contr = chkpt_rd_contr();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  sprim = chkpt_rd_sprim();

  this_atom = snuc[0];
  fprintf(molfile, "  %3d 0\n", this_atom);
  for(i=0; i < nshell; i++) {
    if(this_atom != snuc[i]) {
      this_atom = snuc[i];
      fprintf(molfile, "\n%3d 0\n", this_atom);
    }
    fprintf(molfile, "%s %3d 1.00\n", am_labs[stype[i]-1], snumg[i]);
    for(j=sprim[i]-1; j < sprim[i]-1+snumg[i]; j++) {
      fprintf(molfile, "%20.10f %20.10f\n", exps[j], contr[j]);
    }
  }
  free(snumg);
  free(exps);
  free(contr);
  free(stype);
  free(snuc);

  /** MO's **/
  fprintf(molfile, "[MO]\n");
  fflush(molfile);

  scf_so = chkpt_rd_scf();
  if(scf_so == NULL) {
    printf("You must run cscf before psi2molden.\n");
    exit(PSI_RETURN_SUCCESS);
  }

  nmo = chkpt_rd_nmo();
  nao = chkpt_rd_nao();
  nso = chkpt_rd_nso();
  nirreps = chkpt_rd_nirreps();
  orbspi = chkpt_rd_orbspi();
  clsdpi = chkpt_rd_clsdpi();
  evals = chkpt_rd_evals();
  usotao = chkpt_rd_usotao();
  irr_labs = chkpt_rd_irr_labs();

  /* convert to AO's */
  scf = block_matrix(nao, nmo);
  C_DGEMM('t', 'n', nao, nmo, nso, 1.0, &(usotao[0][0]), nao, &(scf_so[0][0]), nmo, 
	  0.0, &(scf[0][0]), nmo);

  ndocc = 0;
  for(i=0; i < nirreps; i++) ndocc += clsdpi[i];

  /* set up list of MO irreps */
  orbsym = init_int_array(nmo);
  for(h=0,I=0; h < nirreps; h++)
    for(i=0; i < orbspi[h]; i++,I++)
      orbsym[I] = h;

  /* arrange MO's in ascending order */
  mosort(evals, scf, orbsym, nso, nmo);

  for(i=0; i < nmo; i++) {
    fprintf(molfile, " Sym=A\n");
    fprintf(molfile, " Ene=%20.10f\n", evals[i]);
    fprintf(molfile, " Spin=Alpha\n");
    occ = (i < ndocc) ? 2.0 : 0.0;
    fprintf(molfile, " Occup=%3.1f\n", occ);
    for(j=0; j < nso; j++) fprintf(molfile, "%3d %20.12f\n", j+1, scf[j][i]);
  }

  free(orbsym);
  free(orbspi);
  free_block(scf_so);
  free_block(scf);
  free(evals);
  for(i=0; i < nirreps; i++) free(irr_labs[i]);
  free(irr_labs);

  chkpt_close();

  fclose(molfile);

  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

void init_io(int argc, char *argv[])
{
  extern char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0); /* this assumes no cmdline args except filenames */
  ip_cwk_add(":INPUT");
  ip_cwk_add(progid);
  free(progid);
  psio_init(); psio_ipv1_config();
}


void exit_io(void)
{
  psio_done();
  psi_stop();
}

char *gprgid(void)
{
   char *prgid = "PSI2MOLDEN";

   return(prgid);
}
