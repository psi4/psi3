#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <libint.h>
#include <file30.h>
#include <iwl.h>
#include <psifiles.h>
#include <qt.h>

#include <dmalloc.h>

FILE *infile, *outfile;

void init_io(void);
void exit_io(void);
void title(void);

int main(int argc, char *argv[])
{
  int iter, s, t, A, k, l, m, inew, iold, max;
  int i, j, ij, am, atom, shell_length, offset, stat;
  int nirreps, nao, nmo, nso, natom, nshell, noei, nocc, errcod, nfzc;
  int *stype, *snuc, *aostart, *aostop, *ao2atom, *l_length;
  int *clsdpi, *openpi, *orbspi, *dummy, *order, *frdocc;
  double *ss, **S, **scf, **u, **Ctmp, **C, *evals;

  int *orb_order, *orb_boolean;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, Cks, Ckt, **U, **V, **VV, **F;
  double cos4a, alpha, alphamax, alphalast, conv, norm;

  init_io();
  title();

  file30_init();
  nao = file30_rd_nao();
  nmo = file30_rd_nmo();
  nso = file30_rd_nso();
  natom = file30_rd_natom();
  nshell = file30_rd_nshell();
  stype = file30_rd_stype();
  snuc = file30_rd_snuc();
  u = file30_rd_usotao_new();
  nirreps = file30_rd_nirreps();
  clsdpi = file30_rd_clsdpi();
  openpi = file30_rd_openpi();
  orbspi = file30_rd_orbspi();
  scf = file30_rd_scf();
  evals = file30_rd_evals();
  file30_close();

  /* A couple of error traps */
  if(nirreps != 1) {
    fprintf(outfile, "\n\tError: localization is only valid in C1 symmetry!\n");
    exit(2);
  }
  if(nso != nao) {
    fprintf(outfile, "\n\tError: localization is only valid with cartesian polarization!\n");
    exit(2);
  }
  if(openpi[0]) {
    fprintf(outfile, "\n\tError: localization available for closed-shells only!\n");
    exit(2);
  }

  /* Frozen orbital info */
  frdocc = init_int_array(1);
  errcod = ip_int_array("FROZEN_DOCC", frdocc, 1);
  nfzc = frdocc[0];
  free(frdocc);

  /* Compute the length of each AM block */
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=0; l < (LIBINT_MAX_AM); l++) 
    if(l) l_length[l] = l_length[l-1] + l + 1;

  /* Set up the atom->AO and AO->atom lookup arrays */
  aostart = init_int_array(natom);
  aostop = init_int_array(natom);
  for(i=0,atom=-1,offset=0; i < nshell; i++) {
    am = stype[i] - 1;
    shell_length = l_length[am];

    if(atom != snuc[i]-1) {
      if(atom != -1) aostop[atom] = offset-1;
      atom = snuc[i]-1;
      aostart[atom] = offset;
    }

    offset += shell_length;
  }
  aostop[atom] = offset-1;

  ao2atom = init_int_array(nao);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) ao2atom[j] = i;

  /* Get the overlap integrals -- these should be identical to AO S */
  noei = nmo*(nmo+1)/2;
  ss = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_S,ss,noei,0,0,outfile);
  S = block_matrix(nao,nao);
  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++) {
      S[i][j] = S[j][i] = ss[ij];
    }
  free(ss);

  /* transform the MO coefficients to the AO basis */
  C = block_matrix(nao,nmo);
  C_DGEMM('t','n',nao,nmo,nso,1,&(u[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(C[0][0]),nmo);

  /*
  fprintf(outfile, "\tCanonical MO's:\n");
  print_mat(C,nao,nmo,outfile);
  */


  /* Compute nocc --- closed-shells only */
  for(i=0,nocc=0; i < nirreps; i++) nocc += clsdpi[i];

  fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);

  fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
  fprintf(outfile, "\t------------------------------------------------------------\n");

  U = block_matrix(nocc, nocc);
  V = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++) V[i][i] = 1.0;
  VV = block_matrix(nocc, nocc);

  for(iter=0; iter < 100; iter++) {

    P = 0.0;
    for(i=nfzc; i < nocc; i++) {
      for(A=0; A < natom; A++) {
	PiiA = 0.0;

	for(l=aostart[A]; l <= aostop[A]; l++)
	  for(k=0; k < nao; k++) 
	    PiiA += C[k][i] * C[l][i] * S[k][l];

	P += PiiA * PiiA;
      }
    }

    /* Compute 2x2 rotations for Pipek-Mezey localization */
    alphamax = 0.0;
    for(s=nfzc; s < nocc; s++) {
      for(t=nfzc; t < s; t++) {

	Ast = Bst = 0.0;

	for(A=0; A < natom; A++) {

	  Pst = Pss = Ptt = 0.0;
	      
	  for(l=aostart[A]; l <= aostop[A]; l++) {
	    for(k=0; k < nao; k++) {
	      Pst += 0.5 * (C[k][s] * C[l][t] +
			    C[l][s] * C[k][t]) * S[k][l];

	      Pss += C[k][s] * C[l][s] * S[k][l];

	      Ptt += C[k][t] * C[l][t] * S[k][l];
	    }
	  }

	  Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);
	  Bst += Pst * (Pss - Ptt);

	} /* A-loop */

	/* Compute the rotation angle */
	AB = Ast * Ast + Bst * Bst;
	if(fabs(AB) > 0.0) {
	  cos4a = -Ast/sqrt(AB);
	  alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
	}
	else alpha = 0.0;

	/* Keep up with the maximum 2x2 rotation angle */
	alphamax = (alpha > alphamax ? alpha : alphamax);

	Uss = cos(alpha);
	Utt = cos(alpha);
	Ust = sin(alpha);
	Uts = -Ust;

	/* Now do the rotation */
	for(k=0; k < nao; k++) {
	  Cks = C[k][s];
	  Ckt = C[k][t];
	  C[k][s] = Uss * Cks + Ust * Ckt;
	  C[k][t] = Uts * Cks + Utt * Ckt;
	}

	zero_mat(U, nocc, nocc);
	for(i=0; i < nocc; i++) U[i][i] = 1.0;

	U[s][s] = Uss;
	U[t][t] = Utt;
	U[s][t] = Ust;
	U[t][s] = Uts;

	zero_mat(VV, nocc, nocc);
	for(i=0; i < nocc; i++) {
	  for(j=0; j < nocc; j++) {
	    for(k=0; k < nocc; k++) {
	      VV[i][j] += V[i][k] * U[j][k];
	    }
	  }
	}

	for(i=0; i < nocc; i++)
	  for(j=0; j < nocc; j++)
	    V[i][j] = VV[i][j];
	      
      } /* t-loop */
    } /* s-loop */

    conv = fabs(alphamax - alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    if(iter && (conv < 1e-12)) break;
    alphalast = alphamax;

    fflush(outfile);
      
  } /* iter-loop */

  print_mat(V, nocc, nocc, outfile); 

  /* Transform occupied orbital eigenvalues */
  F = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nocc; j++)
      for(k=0; k < nocc; k++) 
	F[i][j] += V[k][i] * evals[k] * V[k][j];

  /*
  fprintf(outfile, "\nTransformed Orbital Energies:\n");
  print_mat(F, nocc, nocc, outfile);
  */

  /* Compute a reordering array based on the diagonal elements of F */
  orb_order = init_int_array(nocc);
  orb_boolean = init_int_array(nocc);
  for(i=0; i < nocc; i++) { orb_order[i] = 0;  orb_boolean[i] = 0; }

  for(i=0,max=0; i < nocc; i++) /* First, find the overall maximum */
    if(fabs(F[i][i]) > fabs(F[max][max])) max = i;

  orb_order[0] = max;  orb_boolean[max] = 1;

  for(i=1; i < nocc; i++) {
    max = 0;
    while(orb_boolean[max]) max++; /* Find an unused max */
    for(j=0; j < nocc; j++) 
      if((fabs(F[j][j]) >= fabs(F[max][max])) && !orb_boolean[j]) max = j;
    orb_order[i] = max; orb_boolean[max] = 1;
  }

  /*
  for(i=0; i < nocc; i++) fprintf(outfile, "%d %d\n", i, orb_order[i]);
  */

  /*
  fprintf(outfile, "\n\tPipek-Mezey Localized MO's (before sort):\n");
  print_mat(C, nao, nmo, outfile);
  */

  /* Now reorder the localized MO's according to F */
  Ctmp = block_matrix(nao,nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nao; j++) Ctmp[j][i] = C[j][i];

  for(i=0; i < nocc; i++) {
    iold = orb_order[i];
    for(j=0; j < nao; j++) C[j][i] = Ctmp[j][iold];
    evals[i] = F[iold][iold];
  }
  free_block(Ctmp);

  /*
  fprintf(outfile, "\n\tPipek-Mezey Localized MO's (after sort):\n");
  print_mat(C, nao, nmo, outfile);
  */

  /* Check MO normalization */
  /*
  for(i=0; i < nmo; i++) {
    norm = 0.0;
    for(j=0; j < nao; j++) 
      for(k=0; k < nao; k++) {
	norm += C[j][i] * C[k][i] * S[j][k];
      }

    fprintf(outfile, "norm[%d] = %20.10f\n", i, norm);
  }
  */

  /* Write the new MO's to file30 */
  file30_init();
  file30_wt_scf(C);
  file30_close();

  free_block(C);
  free(evals);

  fprintf(outfile, "\n\tLocalization of occupied orbitals complete.\n");

  exit_io();
  exit(0);
}

void init_io(void)
{
  char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(progid);

  free(progid);

  psio_init();
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*          LOCAL         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  psio_done();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "LOCAL";

   return(prgid);
}
