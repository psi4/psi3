/*! \defgroup LOCALIZE localize: Localize the orbitals */

/*! 
** \file
** \ingroup LOCALIZE
** \brief Localize the orbitals
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libint/libint.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libqt/qt.h>

extern "C" {
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi { namespace localize {
  void init_io(int argc, char *argv[]);
  void exit_io(void);
  void title(void);
  void do_pipekmezey();
  void do_boysfoster();
}}

int main(int argc, char *argv[])
{
  using namespace psi::localize;
  init_io(argc, argv);
  title();

  char* method;
  int errcod = ip_string("METHOD", &method, 0);
  if (errcod == IPE_OK && strcmp(method,"PIPEK-MEZEY") && strcmp(method,"BOYS-FOSTER") )
    fprintf(outfile, "invalid value for keyword METHOD: must be \"PIPEK-MEZEY\" or \"BOYS-FOSTER\"");
  if (errcod != IPE_OK) // set the default if method is missing
    method = strdup("PIPEK-MEZEY");

  if (strcmp(method,"PIPEK-MEZEY") == 0) {
    fprintf(outfile, "\n\tMethod = Pipek-Mezey\n");
    do_pipekmezey();
  }
  if (strcmp(method,"BOYS-FOSTER") == 0) {
    fprintf(outfile, "\n\tMethod = Boys-Foster\n");
    do_boysfoster();
  }

  free(method);
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

namespace psi { namespace localize {

void do_pipekmezey()
  {
  int iter, s, t, A, k, l, m, p, q, inew, iold, max, max_col, phase_ok, phase_chk;
  int i, j, ij, am, atom, shell_length, offset, stat;
  int nirreps, nao, nmo, nso, natom, nshell, noei, nocc, errcod, nfzc;
  int *stype, *snuc, *aostart, *aostop, *ao2atom, *l_length;
  int *clsdpi, *openpi, *orbspi, *dummy, *order, *frdocc;
  double *ss, **S, **scf, **u, **Ctmp, **C, *evals, **MO_S, **scf_old, **X;

  int *orb_order, *orb_boolean, puream;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, Cks, Ckt, **U, **V, **VV, **F;
  double cos4a, alpha, alphamax, alphalast, conv, norm;
  int print;

  alphalast = 1.0;

  chkpt_init(PSIO_OPEN_OLD);
  puream = chkpt_rd_puream();
  nao = chkpt_rd_nao();
  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();
  natom = chkpt_rd_natom();
  nshell = chkpt_rd_nshell();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  u = chkpt_rd_usotao();
  nirreps = chkpt_rd_nirreps();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  C = chkpt_rd_scf();
  evals = chkpt_rd_evals();
  chkpt_close();

  /* A couple of error traps */
  if(nirreps != 1) {
    fprintf(outfile, "\n\tError: localization is only valid in C1 symmetry!\n");
    exit(PSI_RETURN_FAILURE);
  }
  if(openpi[0]) {
    fprintf(outfile, "\n\tError: localization available for closed-shells only!\n");
    exit(PSI_RETURN_FAILURE);
  }

  /* Frozen orbital info */
  frdocc = get_frzcpi();
  nfzc = frdocc[0];
  free(frdocc);

  /* Compute the length of each AM block */
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=1; l < (LIBINT_MAX_AM); l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

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

  ao2atom = init_int_array(nso);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) ao2atom[j] = i;

  /* Get the overlap integrals -- these should be identical to AO S */
  noei = nso*(nso+1)/2;
  ss = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_S,ss,noei,0,0,outfile);
  S = block_matrix(nso,nso);
  for(i=0,ij=0; i < nso; i++)
    for(j=0; j <= i; j++,ij++) {
      S[i][j] = S[j][i] = ss[ij];
    }
  free(ss);

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
	  for(k=0; k < nso; k++) 
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
	    for(k=0; k < nso; k++) {
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
	alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);

	Uss = cos(alpha);
	Utt = cos(alpha);
	Ust = sin(alpha);
	Uts = -Ust;

	/* Now do the rotation */
	for(k=0; k < nso; k++) {
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

    conv = fabs(alphamax) - fabs(alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    if((iter > 2) && ((fabs(conv) < 1e-12) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);
      
  } /* iter-loop */

  /*  print_mat(V, nocc, nocc, outfile);  */

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
    print_mat(C, nso, nmo, outfile);
  */

  /* Now reorder the localized MO's according to F */
  Ctmp = block_matrix(nso,nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nso; j++) Ctmp[j][i] = C[j][i];

  for(i=0; i < nocc; i++) {
    iold = orb_order[i];
    for(j=0; j < nso; j++) C[j][i] = Ctmp[j][iold];
    evals[i] = F[iold][iold];
  }
  free_block(Ctmp);

  print = 0;
  errcod = ip_boolean("PRINT_MOS", &(print), 0);
  if(print) {
    fprintf(outfile, "\n\tPipek-Mezey Localized MO's (after sort):\n");
    print_mat(C, nso, nmo, outfile);
  }

  /* Check MO normalization */
  /*
    for(i=0; i < nmo; i++) {
    norm = 0.0;
    for(j=0; j < nso; j++) 
    for(k=0; k < nso; k++) {
    norm += C[j][i] * C[k][i] * S[j][k];
    }

    fprintf(outfile, "norm[%d] = %20.10f\n", i, norm);
    }
  */

  /* correct orbital phases for amplitude restarts */
  chkpt_init(PSIO_OPEN_OLD);
  scf_old = chkpt_rd_local_scf();
  chkpt_close();
  if (scf_old != NULL) {
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nmo);
    C_DGEMM('n','n',nso, nmo, nso, 1, &(S[0][0]), nso, &(C[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('t','n',nmo, nmo, nso, 1, &(scf_old[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    /*
    fprintf(outfile, "Approximate Overlap Matrix\n");
    print_mat(MO_S, nmo, nmo, outfile);
    */

    for(p=0; p < nmo; p++) {
      max = 0.0;
      for(q=0; q < nmo; q++) {
	if(fabs(MO_S[p][q]) > max) {
	  max = fabs(MO_S[p][q]); max_col = q;
	}
      }
      if(max_col != p) phase_ok = 0;
    }

    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_phase_check(phase_ok);
    chkpt_close();
    if(phase_ok) {
      for(p=0; p < nmo; p++) {
	if(MO_S[p][p] < 0.0) {
	  for(q=0; q < nso; q++)
	    C[q][p] *= -1.0;
	}
      }
    }

    free_block(MO_S);
    free_block(scf_old);
  }
  free_block(S);

  /* Write the new MO's to chkpt */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_scf(C);
  chkpt_wt_local_scf(C);
  chkpt_close();

  free_block(C);
  free(evals);
}

struct __compare_pairsecond {
    typedef std::pair<int, double> pair_t;
    bool operator()(pair_t A, pair_t B) {
      return A.second < B.second;
    }
};

void do_boysfoster()
{
  chkpt_init(PSIO_OPEN_OLD);
  const int nao = chkpt_rd_nao();
  const int nmo = chkpt_rd_nmo();
  const int nso = chkpt_rd_nso();
  const int nirreps = chkpt_rd_nirreps();
  int* clsdpi = chkpt_rd_clsdpi();
  int* openpi = chkpt_rd_openpi();
  int* orbspi = chkpt_rd_orbspi();
  double** C_so = chkpt_rd_scf();
  double** usotao = chkpt_rd_usotao();
  double** C_ao = block_matrix(nao, nmo);
  double* evals = chkpt_rd_evals();
  mmult(usotao, 1, C_so, 0, C_ao, 0, nao, nso, nmo, 0);
  chkpt_close();

  /* A couple of error traps */
  if(nirreps != 1) {
    fprintf(outfile, "\n\tError: localization is only valid in C1 symmetry!\n");
    exit(PSI_RETURN_FAILURE);
  }
  if(openpi[0]) {
    fprintf(outfile, "\n\tError: localization available for closed-shells only!\n");
    exit(PSI_RETURN_FAILURE);
  }


  /* Get the one-electron integrals */
  const int ntri_ao = nao*(nao+1)/2;
  double* ss = init_array(ntri_ao);
  double* mmx = init_array(ntri_ao);
  double* mmy = init_array(ntri_ao);
  double* mmz = init_array(ntri_ao);
  int stat;
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_S,ss,ntri_ao,0,0,outfile);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MX,mmx,ntri_ao,0,0,outfile);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MY,mmy,ntri_ao,0,0,outfile);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MZ,mmz,ntri_ao,0,0,outfile);
  double** S_ao = block_matrix(nao,nao);
  double** MX_ao = block_matrix(nao,nao);
  double** MY_ao = block_matrix(nao,nao);
  double** MZ_ao = block_matrix(nao,nao);
  int ij = 0;
  for(int i=0; i < nao; i++) {
    for(int j=0; j <= i; j++,ij++) {
      S_ao[i][j] = S_ao[j][i] = ss[ij];
      MX_ao[i][j] = MX_ao[j][i] = mmx[ij];
      MY_ao[i][j] = MY_ao[j][i] = mmy[ij];
      MZ_ao[i][j] = MZ_ao[j][i] = mmz[ij];
    }
  }
  free(ss);
  free(mmx);
  free(mmy);
  free(mmz);
  // transform to occupied MO basis
  int nocc = 0;
  for(int i=0; i < nirreps; i++) nocc += clsdpi[i];
  double** Co = block_matrix(nao, nocc);
  // extract occupied orbital coefficients and transform to AO basis
  {
    int oo = 0;
    int ii = 0;
    for(int h=0; h < nirreps; h++) {
      for(int o=0; o<clsdpi[h]; ++o, ++oo, ++ii) {
        for(int ao=0; ao<nao; ++ao)
          Co[ao][oo] = C_ao[ao][ii];
      }
      ii += orbspi[h] - clsdpi[h];
    }
  }
  double** MX = block_matrix(nocc, nocc);
  double** MY = block_matrix(nocc, nocc);
  double** MZ = block_matrix(nocc, nocc);
  {
    double** tmp = block_matrix(nocc, nao);
    mmult(Co, 1, MX_ao, 0, tmp, 0, nocc, nao, nao, 0);
    mmult(tmp, 0, Co, 0, MX, 0, nocc, nao, nocc, 0);
    mmult(Co, 1, MY_ao, 0, tmp, 0, nocc, nao, nao, 0);
    mmult(tmp, 0, Co, 0, MY, 0, nocc, nao, nocc, 0);
    mmult(Co, 1, MZ_ao, 0, tmp, 0, nocc, nao, nao, 0);
    mmult(tmp, 0, Co, 0, MZ, 0, nocc, nao, nocc, 0);
    free_block(tmp);
  }
  fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n", nocc);

  // compute D = \sum_i <i|r|i>^2
  double D = 0.0;
  for (int i = 0; i < nocc; ++i) {
    D += MX[i][i]*MX[i][i] + MY[i][i]*MY[i][i] + MZ[i][i]*MZ[i][i];
  }
  fprintf(outfile, "\tinitial ext. localization: %lf\n\n", D);

  fprintf(outfile, "\tIter     Ext. Localization   Max. Rotation Angle       Conv\n");
  fprintf(outfile, "\t------------------------------------------------------------\n");

  // U will keep track of rotations
  // C(local) = C(canon) * U
  double** U = block_matrix(nocc, nocc);
  for(int i=0; i < nocc; i++) U[i][i] = 1.0;

  double alphalast = 0.0;
  const unsigned int max_niter = 100;
  for(int iter=0; iter < max_niter; iter++) {

    double alphamax = 0.0;

    /* Compute 2x2 rotations that minimize orbital extent
       (or, equivalently, maximize sum of extents from the origin) */
    for(int s=0; s < nocc; s++) {
      for(int t=0; t < s; t++) {

        const double Ast = MX[s][t]*MX[s][t] + MY[s][t]*MY[s][t] + MZ[s][t]*MZ[s][t] -
                           0.25 * ( (MX[s][s] - MX[t][t])*(MX[s][s] - MX[t][t]) +
                                    (MY[s][s] - MY[t][t])*(MY[s][s] - MY[t][t]) +
                                    (MZ[s][s] - MZ[t][t])*(MZ[s][s] - MZ[t][t])
                                  );

        const double Bst = (MX[s][s] - MX[t][t])*MX[s][t] +
                           (MY[s][s] - MY[t][t])*MY[s][t] +
                           (MZ[s][s] - MZ[t][t])*MZ[s][t];

        /* Compute the rotation angle */
        const double ABst = Ast * Ast + Bst * Bst;
        double alpha = 0.0;
        if(fabs(ABst) > 0.0) {
          const double cos4a = -Ast/sqrt(ABst);
          alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
        }

        /* Keep up with the maximum 2x2 rotation angle */
        alphamax = std::max(alpha,alphamax);

        const double Uss = cos(alpha);
        const double Utt = cos(alpha);
        const double Uts = sin(alpha);
        const double Ust = -Uts;

        /* Rotate orbital cofficients */
        for(int k=0; k < nocc; k++) {
          const double Uks = U[k][s];
          const double Ukt = U[k][t];
          U[k][s] = Uks * Uss + Ukt * Uts;
          U[k][t] = Uks * Ust + Ukt * Utt;
        }
        /* Rotate dipole moment matrices */
        for(int xyz = 0; xyz<3; ++xyz) {
          double** M;
          switch(xyz) {
            case 0: M = MX; break;
            case 1: M = MY; break;
            case 2: M = MZ; break;
          }

          for(int k=0; k < nocc; k++) {
            const double Mks = M[k][s];
            const double Mkt = M[k][t];
            M[k][s] = Mks * Uss + Mkt * Uts;
            M[k][t] = Mks * Ust + Mkt * Utt;
          }
          for(int k=0; k < nocc; k++) {
            const double Msk = M[s][k];
            const double Mtk = M[t][k];
            M[s][k] = Uss * Msk + Uts * Mtk;
            M[t][k] = Ust * Msk + Utt * Mtk;
          }
        }

      } /* t-loop */
    } /* s-loop */

    // compute D = \sum_i <i|r|i>^2
    double D = 0.0;
    for (int i = 0; i < nocc; ++i) {
      D += MX[i][i]*MX[i][i] + MY[i][i]*MY[i][i] + MZ[i][i]*MZ[i][i];
    }

    const double conv = fabs(alphamax) - fabs(alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, D, alphamax, conv);
    if((iter > 2) && ((fabs(conv) < 1e-12) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);

  } /* iter-loop */

#if 0
  fprintf(outfile, "Localization transformation (U):\n");
  print_mat(U, nocc, nocc, outfile);
#endif

  /* Transform occupied orbital eigenvalues */
  double** F = block_matrix(nocc, nocc);
  for(int i=0; i < nocc; i++)
    for(int j=0; j < nocc; j++)
      for(int k=0; k < nocc; k++)
        F[i][j] += U[k][i] * evals[k] * U[k][j];

  fprintf(outfile, "\nFock matrix in the localized MO basis:\n");
  print_mat(F, nocc, nocc, outfile);

  /* Compute a reordering array based on the diagonal elements of F */
  std::vector< std::pair<int, double> >  F_diag(nocc);
  for(int i=0; i<nocc; ++i) F_diag[i] = std::make_pair(i, F[i][i]);
  __compare_pairsecond comparer;
  std::stable_sort(F_diag.begin(), F_diag.end(), comparer);

  /* Now reorder the localized MO's according to F */
  double** U_reord = block_matrix(nocc, nocc);
  for(int i=0; i < nocc; i++) {
    const int iold = F_diag[i].first;
    for(int j=0; j<nocc; ++j)
      U_reord[j][i] = U[j][iold];
  }
  double** C_ao_local = block_matrix(nao, nmo);  std::copy(C_ao[0], C_ao[0]+nao*nmo, C_ao_local[0]);
  for(int i=0; i<nocc; ++i) {
    for(int j=0; j<nao; ++j) {
      double c_ji = 0.0;
      for(int k=0; k<nocc; ++k)
        c_ji += C_ao[j][k] * U_reord[k][i];
      C_ao_local[j][i] = c_ji;
    }
  }
  double** C_so_local = block_matrix(nso, nmo);  std::copy(C_so[0], C_so[0]+nso*nmo, C_so_local[0]);
  for(int i=0; i<nocc; ++i) {
    for(int j=0; j<nso; ++j) {
      double c_ji = 0.0;
      for(int k=0; k<nocc; ++k)
        c_ji += C_so[j][k] * U_reord[k][i];
      C_so_local[j][i] = c_ji;
    }
  }

  int print = 0;
  int errcod = ip_boolean("PRINT_MOS", &(print), 0);
  if(print) {
    fprintf(outfile, "\n\tBoys-Foster Canonical MO's:\n");
    print_mat(C_ao, nao, nmo, outfile);
    fprintf(outfile, "\n\tBoys-Foster Localized MO's:\n");
    print_mat(C_ao_local, nao, nmo, outfile);

    /* Check MO orthonormality */
    double** tmp = block_matrix(nao, nmo);
    mmult(S_ao, 0, C_ao_local, 0, tmp, 0, nao, nao, nmo, 0);
    double** S_mo = block_matrix(nmo, nmo);
    mmult(C_ao_local, 1, tmp, 0, S_mo, 0, nmo, nao, nmo, 0);
    fprintf(outfile, "\n\tOverlap matrix in localized MO basis (should be unit matrix)");
    print_mat(S_mo, nmo, nmo, outfile);

    mmult(C_ao, 1, tmp, 0, S_mo, 0, nmo, nao, nmo, 0);
    fprintf(outfile, "\n\tOverlap matrix between localized and canonical MOs");
    print_mat(S_mo, nmo, nmo, outfile);

    free_block(tmp);
    free_block(S_mo);
  }

  /* correct orbital phases for amplitude restarts */
  chkpt_init(PSIO_OPEN_OLD);
  double** C_so_local_chkpt = chkpt_rd_local_scf();
  chkpt_close();
  if (C_so_local_chkpt != NULL) {
    double** C_ao_local_chkpt = block_matrix(nao, nmo);
    mmult(usotao, 1, C_so_local_chkpt, 0, C_ao_local_chkpt, 0, nao, nso, nmo, 0);

    /* compute overlap between chkpt and computed localized orbitals -- overlap matrix should be "diagonal" */
    double** tmp = block_matrix(nao, nmo);
    mmult(S_ao, 0, C_ao_local, 0, tmp, 0, nao, nao, nmo, 0);
    double** S_mo = block_matrix(nmo, nmo);
    mmult(C_ao_local_chkpt, 1, tmp, 0, S_mo, 0, nmo, nao, nmo, 0);

    int phase_ok = 1;
    for(int p=0; p < nmo; p++) {
      double max = 0.0;
      int max_col;
      for(int q=0; q < nmo; q++) {
        if(fabs(S_mo[p][q]) > max) {
          max = fabs(S_mo[p][q]); max_col = q;
        }
      }
      if(max_col != p) phase_ok = 0;
    }

    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_phase_check(phase_ok);
    chkpt_close();
    if(phase_ok) {
      for(int p=0; p < nmo; p++) {
        if(S_mo[p][p] < 0.0) {
          for(int q=0; q < nso; q++)
            C_so_local[q][p] *= -1.0;
          for(int q=0; q < nao; q++)
            C_ao_local[q][p] *= -1.0;
        }
      }
    }

    free_block(tmp);
    free_block(S_mo);
    free_block(C_so_local_chkpt);
    free_block(C_ao_local_chkpt);
  }

  /* Write the new MO's to chkpt */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_scf(C_so_local);
  chkpt_wt_local_scf(C_so_local);
  chkpt_close();

}

extern "C" {
extern const char *gprgid();
}

void init_io(int argc, char * argv[])
{
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

  psio_init(); psio_ipv1_config();
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
  tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}

}} // namespace psi::localize

extern "C" const char *gprgid()
{
   const char *prgid = "LOCAL";

   return(prgid);
}
