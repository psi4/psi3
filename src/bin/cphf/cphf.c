/*
** CPHF: Program to solve the Coupled Perturbed Hartree-Fock Equations
**
** Be careful!  We assume RHF cases ONLY for now!!!!
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <psio.h>
#include <iwl.h>
#include <file30.h>
#include <qt.h>
#include <psifiles.h>
#include <masses.h>
#include <physconst.h>

FILE *infile, *outfile;
int *ioff;
#define IOFF_MAX 32641

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void init_io(void);
void exit_io(void);
void title(void);
void init_ioff(void);

int main(int argc, char *argv[])
{
  int h, p, q, k, l, i, j, a, b, m;
  int pq, kl, pqkl, pk, ql, pkql, pl, qk, plqk;
  int pj, jq, jl, pjql, pqjl, pljq;
  int r, s, rs, pqrs;
  int psym, qsym, ksym, lsym, rsym, ssym, pqsym;
  int ai, bj, ab, ij, aj, ib, aibj, abij, ajib;
  int AI, BJ, IJ, BK;
  int bk, jk, ik, jb, ijbk, ibjk, ikjb;
  int asym, bsym, isym, jsym;
  int pfirst, qfirst, kfirst, lfirst, afirst, bfirst, ifirst, jfirst, rfirst, sfirst;
  int plast, qlast, klast, llast, alast, blast, ilast, jlast,rlast, slast;
  int natom, nmo, nso, ntri, ntei, nirreps, num_ai, num_ij;
  int *orbspi, *clsdpi, *uoccpi, ndocc, nuocc;
  int *first, *last, *ofirst, *olast, *vfirst, *vlast;
  int foffset, loffset;
  int coord, coord_a, coord_b;
  int error, *ipiv;
  char *label;
  double *inbuf, *evals, *ints, value, det;
  double ***F, ***S, ***B, **A, **AA, **B0, ***U;
  double **hessian, *eta, **Acopy;
  double **M, **evecs, *km, k_convert, cm_convert;
  double *zvals, **temp;
  FILE *file15;

  init_io();
  title();
  init_ioff();

  timer_init();
  timer_on("CPHF Main");

  file30_init();
  natom = file30_rd_natom();
  nmo = file30_rd_nmo();
  nso = file30_rd_nso();
  evals = file30_rd_evals();
  nirreps = file30_rd_nirreps();
  orbspi = file30_rd_orbspi();
  clsdpi = file30_rd_clsdpi();
  zvals = file30_rd_zvals();
  file30_close();

  /*
  fprintf(outfile, "Orbital Eigenvalues:\n");
  for(i=0; i < nmo; i++) fprintf(outfile, "%d %20.12f\n", i, evals[i]);
  */

  ntri = nmo * (nmo + 1)/2;

  uoccpi = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) uoccpi[h] = orbspi[h] - clsdpi[h];

  ndocc = 0; nuocc = 0;
  for(h=0; h < nirreps; h++) { ndocc += clsdpi[h]; nuocc += uoccpi[h]; }

  num_ai = ndocc * nuocc; 
  num_ij = ndocc * (ndocc + 1)/2;

  /* Build the first and last lookup arrays */
  first = init_int_array(nirreps);
  last = init_int_array(nirreps);
  foffset = 0;
  loffset = orbspi[0] - 1;
  first[0] = foffset;
  last[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += orbspi[h-1];
    loffset += orbspi[h];
    first[h] = foffset;
    last[h] = loffset;
  }

  ofirst = init_int_array(nirreps);
  olast = init_int_array(nirreps);
  foffset = 0;
  loffset = clsdpi[0] - 1;
  ofirst[0] = foffset;
  olast[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += orbspi[h-1];
    loffset += uoccpi[h-1] + clsdpi[h];
    ofirst[h] = foffset;
    olast[h] = loffset;
  }

  vfirst = init_int_array(nirreps);
  vlast = init_int_array(nirreps);
  foffset = clsdpi[0];
  loffset = orbspi[0] - 1;
  vfirst[0] = foffset;
  vlast[0] = loffset;
  for(h=1; h < nirreps; h++) {
    foffset += uoccpi[h-1] + clsdpi[h];
    loffset += orbspi[h];
    vfirst[h] = foffset;
    vlast[h] = loffset;
  }

  /* Allocate space for the F, S, and B vectors */
  F = (double ***) malloc(natom*3 * sizeof(double **));
  S = (double ***) malloc(natom*3 * sizeof(double **));
  B = (double ***) malloc(natom*3 * sizeof(double **));
  for(coord=0; coord < natom*3; coord++) {
    F[coord] = block_matrix(nmo, nmo);
    S[coord] = block_matrix(nmo, nmo);
    B[coord] = block_matrix(nmo, nmo);
  }

  /* Grab the MO-basis overlap and Fock derivative integrals from disk */
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  inbuf = init_array(ntri);
  for(coord=0; coord < natom*3; coord++) {

    sprintf(label, "MO-basis Fock Derivs (%d)", coord);
    iwl_rdone(PSIF_OEI, label, inbuf, ntri, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0, ij=0; i < nmo; i++)
      for(j=0; j <= i; j++, ij++)
	F[coord][i][j] = F[coord][j][i] = inbuf[ij];

/*
    fprintf(outfile, "F[%d] Deriv (MO):\n", coord);
    print_mat(F[coord], nmo, nmo, outfile);
*/
  }

  for(coord=0; coord < natom*3; coord++) {
    sprintf(label, "MO-basis Overlap Derivs (%d)", coord);
    iwl_rdone(PSIF_OEI, label, inbuf, ntri, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0, ij=0; i < nmo; i++)
      for(j=0; j <= i; j++, ij++)
	S[coord][i][j] = S[coord][j][i] = inbuf[ij];

/*
    fprintf(outfile, "S[%d] Deriv (MO):\n", coord);
    print_mat(S[coord], nmo, nmo, outfile);
*/
  }
  free(inbuf);
  free(label);

  /***** Build the B vectors *****/

  /* Fock and overlap derivative components */
  for(coord=0; coord < natom*3; coord++) {

    for(i=0; i < nmo; i++) {
      for(j=0; j < nmo; j++) {
	B[coord][i][j] = F[coord][i][j] - S[coord][i][j] * evals[j];
      }
    }
  }

  /*
  for(coord=0; coord < natom*3; coord++) {
    fprintf(outfile, "\nB0[%d] Vector (MO) (Fa & Sa only) :\n", coord);
    print_mat(B[coord], nmo, nmo, outfile);
  }
  */

  /* Grab the two-electron integrals */
  ntei = ntri * (ntri + 1)/2;
  ints = init_array(ntei);
  iwl_rdtwo(PSIF_MO_TEI, ints, ioff, nmo, 0, 0, 0, outfile);

  /* Two-electron integral components of B vectors */
  for(psym=0; psym < nirreps; psym++) {
    pfirst = first[psym];
    plast = last[psym];

    for(p=pfirst; p <= plast; p++) {

      for(qsym=0; qsym < nirreps; qsym++) {
	qfirst = first[qsym];
	qlast = last[qsym];

	for(q=qfirst; q <= qlast; q++) {
	  pq = INDEX(p,q);

	  for(ksym=0; ksym < nirreps; ksym++) {
	    kfirst = ofirst[ksym];
	    klast = olast[ksym];

	    for(k=kfirst; k <= klast; k++) {
	      pk = INDEX(p,k);
	      qk = INDEX(q,k);

	      for(lsym=0; lsym < nirreps; lsym++) {

		lfirst = ofirst[lsym];
		llast = olast[lsym];

		for(l=lfirst; l <= llast; l++) {
		  kl = INDEX(k,l);
		  ql = INDEX(q,l);
		  pl = INDEX(p,l);

		  pqkl = INDEX(pq,kl);
		  pkql = INDEX(pk,ql);
		  plqk = INDEX(pl,qk);

		  value = 2.0 * ints[pqkl] - ints[pkql];

		  for(coord=0; coord < natom*3; coord++) {
		    B[coord][p][q] -= S[coord][k][l] * value;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /*
  for(coord=0; coord < natom*3; coord++) {
    fprintf(outfile, "\nB0[%d] Matrix (MO):\n", coord);
    print_mat(B[coord], nmo, nmo, outfile);
  }
  */

  /* Sort the B vectors */
  B0 = block_matrix(natom*3, num_ai);
  for(coord=0; coord < natom*3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym = 0; isym < nirreps; isym++) {
	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    B0[coord][AI] = B[coord][a][i];
	  }
	}
      }
    }
  }

  /* Build the MO Hessian */
  A = block_matrix(num_ai, num_ai);
  for(asym=0,AI=0; asym < nirreps; asym++) {

    afirst = vfirst[asym];
    alast = vlast[asym];

    for(a=afirst; a <= alast; a++) {

      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];

	for(i=ifirst; i <= ilast; i++,AI++) {
	  ai = INDEX(a,i);

	  for(bsym=0,BJ=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
	      ab = INDEX(a,b);
	      ib = INDEX(i,b);

	      for(jsym=0; jsym < nirreps; jsym++) {

		jfirst = ofirst[jsym];
		jlast = olast[jsym];

		for(j=jfirst; j <= jlast; j++,BJ++) {
		  bj = INDEX(b,j);
		  ij = INDEX(i,j);
		  aj = INDEX(a,j);

		  aibj = INDEX(ai,bj);
		  abij = INDEX(ab,ij);
		  ajib = INDEX(aj,ib);

		  A[AI][BJ] = (a==b) * (i==j) * (evals[i] - evals[a]);
		  A[AI][BJ] -= 4.0 * ints[aibj] - ints[abij] - ints[ajib];

		}
	      }
	    }
	  }
	}
      }
    }
  }

  /* Make a copy of A */
  Acopy = block_matrix(num_ai, num_ai);
  for(AI=0; AI < num_ai; AI++)
    for(BJ=0; BJ < num_ai; BJ++)
      Acopy[AI][BJ] = A[AI][BJ];


/*
  fprintf(outfile, "\nMO Hessian:\n");
  print_mat(A, num_ai, num_ai, outfile);
*/

  /* Build the Mixed Hessian (dependent x independent) */
  AA = block_matrix(num_ij, num_ai);
  for(isym=0,IJ=0; isym < nirreps; isym++) {

    ifirst = ofirst[isym];
    ilast = olast[isym];

    for(i=ifirst; i <= ilast; i++) {

      for(jsym=0; jsym < nirreps; jsym++) {

	jfirst = ofirst[jsym];
	jlast = olast[jsym];

	for(j=jfirst; (j <= jlast) && (j <= i); j++,IJ++) {
	  ij = INDEX(i,j);

	  if(i==j) continue;

	  for(bsym=0,BK=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
	      ib = INDEX(i,b);
	      jb = INDEX(j,b);

	      for(ksym=0; ksym < nirreps; ksym++) {

		kfirst = ofirst[ksym];
		klast = olast[ksym];

		for(k=kfirst; k <= klast; k++,BK++) {
		  bk = INDEX(b,k);
		  ik = INDEX(i,k);
		  jk = INDEX(j,k);

		  ijbk = INDEX(ij,bk);
		  ibjk = INDEX(ib,jk);
		  ikjb = INDEX(ik,jb);

		  AA[IJ][BK] = 4.0 * ints[ijbk] - ints[ibjk] - ints[ikjb];

		}
	      }
	    }
	  }
	}

      }
    }
  }

  /*
  fprintf(outfile, "\nMixed MO Hessian:\n");
  print_mat(AA, num_ij, num_ai, outfile);
  */

  ipiv = init_int_array(num_ai);

  /* Solve the CPHF equations */
  for(coord=0; coord < natom*3; coord++) {

/*
    fprintf(outfile, "\nB0[%d] vector before flin:\n", coord);
    for(ai=0; ai < num_ai; ai++) fprintf(outfile, "%d %20.12f\n", ai, B0[coord][ai]);
*/

    /* flin(A, B0[coord], num_ai, 1, &det); */
    /* pople(A, B0[coord], num_ai, 1, 1e-10, outfile, 0); */
    error = C_DGESV(num_ai, 1, &(A[0][0]), num_ai, &(ipiv[0]), &(B0[coord][0]), num_ai);
 /*   fprintf(outfile, "\nerror = %d\n", error); */

    /*    fprintf(outfile, "\nDeterminant: %20.12f\n", det); */

/*
    fprintf(outfile, "\nU[%d] vector after flin:\n", coord);
    for(ai=0; ai < num_ai; ai++) fprintf(outfile, "%d %20.12f\n", ai, B0[coord][ai]);
*/

    /* Recopy A */
    for(AI=0; AI < num_ai; AI++)
      for(BJ=0; BJ < num_ai; BJ++)
	A[AI][BJ] = Acopy[AI][BJ];
  }

  /* Sort the U matrices to full MO form */
  U = (double ***) malloc(natom*3 * sizeof(double **));
  for(coord=0; coord < natom*3; coord++) U[coord] = block_matrix(nmo, nmo);
  for(coord=0; coord < natom*3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym=0; isym < nirreps; isym++) {

	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    U[coord][a][i] = B0[coord][AI];
	  }
	}
      }
    }
  }

  /* Add the dependent pairs on the lower triangle */
  for(isym=0,IJ=0; isym < nirreps; isym++) {
    ifirst = ofirst[isym];
    ilast = olast[isym];

    for(i=ifirst; i <= ilast; i++) {

      for(jsym=0; jsym < nirreps; jsym++) {
	jfirst = ofirst[jsym];
	jlast = olast[jsym];

	for(j=jfirst; (j <= jlast) && (j <= i); j++,IJ++) {

	  if(i==j) continue;

	  for(coord=0; coord < natom*3; coord++) U[coord][i][j] = -0.5 * S[coord][i][j];

	  /*

	  for(coord=0; coord < natom*3; coord++) U[coord][i][j] = B[coord][i][j];

	  for(bsym=0,BK=0; bsym < nirreps; bsym++) {
	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {

	      for(ksym=0; ksym < nirreps; ksym++) {
		kfirst = ofirst[ksym];
		klast = olast[ksym];

		for(k=kfirst; k <= klast; k++,BK++) {

		  for(coord=0; coord < natom*3; coord++) 
		    U[coord][i][j] += AA[IJ][BK] * U[coord][b][k];
		}
	      }
	    }
	  }
	  */

	  /*
	  for(coord=0; coord < natom*3; coord++) U[coord][i][j] *= 1.0/(evals[j] - evals[i]);
	  */
	}
      }
    }
  }


  /* Add the upper triangle --- careful here! */
  for(coord=0; coord < natom*3; coord++) {

    /* dependent pairs */
    for(isym=0; isym < nirreps; isym++) {
      ifirst = ofirst[isym];
      ilast = olast[isym];
      for(i=ifirst; i <= ilast; i++) {

	U[coord][i][i] -= 0.5 * S[coord][i][i];

	for(jsym=0; jsym < nirreps; jsym++) {
	  jfirst = ofirst[jsym];
	  jlast = olast[jsym];

	  for(j=jfirst; (j <= jlast) && (j < i); j++)
	    U[coord][j][i] = -U[coord][i][j] - S[coord][i][j];

	}
      }
    }

    /* independent pairs */
    for(asym=0; asym < nirreps; asym++) {
      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym=0; isym < nirreps; isym++) {
	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++)
	    U[coord][i][a] = -U[coord][a][i] - S[coord][a][i];

	}
      }
    }
  }

  /* Print the Ua matrices */
/*
  for(coord=0; coord < natom*3; coord++) {
    fprintf(outfile, "\nU[%d] Matrix (MO):\n", coord);
    print_mat(U[coord], nmo, nmo, outfile);
  }
*/

  /*** Compute the SCF second derivative ***/
  hessian = block_matrix(natom*3, natom*3);

  /* grab skeleton derivatives from cints --deriv2 **/
  psio_open(PSIF_DERINFO, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_DERINFO, "Skeleton Hessian", (char *) hessian[0], natom*3*natom*3*sizeof(double));
  psio_close(PSIF_DERINFO, 1);

/*
  fprintf(outfile, "Skeleton Contribution to SCF Molecular Hessian:\n");
  print_mat(hessian, natom*3, natom*3, outfile);
*/

  eta = init_array(nmo);
  for(coord_a=0; coord_a < natom*3; coord_a++) {
    for(coord_b=0; coord_b < natom*3; coord_b++) {

      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];
	for(i=ifirst; i <= ilast; i++) {
	  eta[i] = 0.0;
	  for(m=0; m < nmo; m++) {
	    eta[i] += (2.0 * U[coord_a][i][m] * U[coord_b][i][m] - 
		       2.0 * S[coord_a][i][m] * S[coord_b][i][m]);
	  }
	}
      }

      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];
	for(i=ifirst; i <= ilast; i++) 
	  hessian[coord_a][coord_b] -= 2.0 * eta[i] * evals[i];
      }

      for(psym=0; psym < nirreps; psym++) {

	pfirst = first[psym];
	plast = last[psym];

	for(p=pfirst; p <= plast; p++) {

	  for(jsym=0; jsym < nirreps; jsym++) {
	    jfirst = ofirst[jsym];
	    jlast = olast[jsym];

	    for(j=jfirst; j <= jlast; j++) {

	      pj = INDEX(p,j);

	      hessian[coord_a][coord_b] += 4.0 * (U[coord_b][p][j] * F[coord_a][p][j] + 
						  U[coord_a][p][j] * F[coord_b][p][j]);

	      hessian[coord_a][coord_b] += 4.0 * U[coord_a][p][j] * U[coord_b][p][j] * evals[p];

	      for(qsym=0; qsym < nirreps; qsym++) {
		qfirst = first[qsym];
		qlast = last[qsym];

		for(q=qfirst; q <= qlast; q++) {

		  pq = INDEX(p,q);
		  jq = INDEX(j,q);

		  for(lsym=0; lsym < nirreps; lsym++) {
		    lfirst = ofirst[lsym];
		    llast = olast[lsym];

		    for(l=lfirst; l <= llast; l++) {

		      ql = INDEX(q,l);
		      jl = INDEX(j,l);
		      pl = INDEX(p,l);

		      pjql = INDEX(pj,ql);
		      pqjl = INDEX(pq,jl);
		      pljq = INDEX(pl,jq);

		      hessian[coord_a][coord_b] += 4.0 * U[coord_a][p][j] * U[coord_b][q][l] *
			(4.0 * ints[pjql] - ints[pqjl] - ints[pljq]);
		    }
		  }
		}
	      }

	    }
	  }
	}
      }

    }
  }
  free(eta);

  fprintf(outfile, "SCF Molecular Hessian:\n");
  print_mat(hessian, natom*3, natom*3, outfile);

  /* write the hessian to file15 in the PSI2 standard format */
  file15 = fopen("file15.dat", "w");
  fprintf(file15, "%5d%5d\n", natom, natom*6);
  for(i=0; i < natom*3; i++) {
    for(j=0; j < natom; j++) {
      fprintf(file15, "%20.10f%20.10f%20.10f\n", 
	      hessian[i][j*3], hessian[i][j*3+1], hessian[i][j*3+2]);
    }
  }
  fclose(file15);

  /**** compute the vibrational frequencies ****/

  /* mass-weight the hessian */
  M = block_matrix(natom*3, natom*3);
  for(i=0; i < natom; i++) {
    for(j=0; j < 3; j++)  {
      M[i*3+j][i*3+j] = 1/sqrt(an2masses[(int) zvals[i]]);
    }
  }

  temp = block_matrix(natom*3,natom*3);
  C_DGEMM('n','n', natom*3, natom*3, natom*3, 1.0, &(M[0][0]), natom*3,
          &(hessian[0][0]), natom*3, 0.0, &(temp[0][0]), natom*3);
  C_DGEMM('n','n', natom*3, natom*3, natom*3, 1.0, &(temp[0][0]), natom*3,
	  &(M[0][0]), natom*3, 0.0, &(hessian[0][0]), natom*3);
  free_block(temp);
  free_block(M);

  /*

  fprintf(outfile, "Mass-Weighting matrix:\n");
  print_mat(M, natom*3, natom*3, outfile);

  fprintf(outfile, "Mass-Weighted Hessian:\n");
  print_mat(hessian, natom*3, natom*3, outfile);
  */

  /* diagonalize mass-weighted hessian */
  km = init_array(natom*3);
  evecs = block_matrix(natom*3, natom*3);
  sq_rsp(natom*3, natom*3, hessian, km, 1, evecs, 1e-12);

  /* compute the frequencies and spit them out in a nice table */
  fprintf(outfile, "\n\tHarmonic Vibrational Frequencies (cm-1):\n");
  fprintf(outfile,   "\t----------------------------------------\n");
  k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for(i=natom*3-1; i >= 0; i--) {
    if(km[i] < 0.0)
      fprintf(outfile, "\t  %3d   %17.3fi\n", i, cm_convert * sqrt(-k_convert * km[i]));
    else
      fprintf(outfile, "\t  %3d   %17.3f\n", i, cm_convert * sqrt(k_convert * km[i]));
  }

  for(coord=0; coord < natom*3; coord++) {
    free_block(F[coord]);
    free_block(S[coord]);
    free_block(B[coord]);
    free_block(U[coord]);
  }
  free(F); free(S); free(B); free(U);

  free_block(B0); free_block(A); free_block(AA);

  free(first); free(last);
  free(ofirst); free(olast);
  free(vfirst); free(vlast);

  free_block(hessian);
  free_block(evecs);
  free(km);
  free(zvals);

  free(ints);
  free(evals);
  free(orbspi);
  free(clsdpi);
  free(uoccpi);
  free(ioff);

  timer_off("CPHF Main");
  timer_done();

  exit_io();
  exit(0);
}

void init_io(void)
{
  int i;
  char *gprgid(void);
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
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*          CPHF          *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  psio_done();

  free_ptrs();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid(void)
{
   char *prgid = "CPHF";

   return(prgid);
}

void init_ioff(void)
{
  int i;
  ioff = (int *) malloc(IOFF_MAX * sizeof(int));
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) {
      ioff[i] = ioff[i-1] + i;
    }
}
