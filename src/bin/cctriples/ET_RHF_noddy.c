#include <stdio.h>
#include <math.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double ***init_3d_array(int p, int q, int r);
void free_3d_array(double ***A, int p, int q);

double ET_RHF_noddy(void)
{
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gij, Gji, Gjk, Gkj, Gik, Gki;
  int Gab, Gba, Gac, Gca, Gbc, Gcb;
  int I, J, K, A, B, C, D, L;
  int i, j, k, a, b, c, d, l;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int id, kd, jd, ad, bd, cd;
  int la, lb, lc, il, jl, kl;
  int nv;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double value, dijk, value1, value2, denom, ET;
  double ***W, ***V, ***F;
  double t_kjcd, t_jkbd, t_jibd, t_ijad, t_ikad, t_kicd;
  double F_idab, F_idac, F_kdca, F_kdcb, F_jdbc, F_jdba;
  double t_ilab, t_ilac, t_klca, t_klcb, t_jlbc, t_jlba;
  double E_jklc, E_kjlb, E_ijlb, E_jila, E_kila, E_iklc;
  double t_ia, t_jb, t_kc, D_jkbc, D_ikac, D_ijab;
  dpdbuf4 T2, Fints, Eints, Dints;
  dpdfile2 fIJ, fAB, T1;
  double chksum;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);

  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2, h);
    dpd_buf4_mat_irrep_rd(&T2, h);

    dpd_buf4_mat_irrep_init(&Fints, h);
    dpd_buf4_mat_irrep_rd(&Fints, h);

    dpd_buf4_mat_irrep_init(&Eints, h);
    dpd_buf4_mat_irrep_rd(&Eints, h);

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);
  }

  for(h=0,nv=0; h < nirreps; h++) nv += virtpi[h];

  ET = 0.0;

  F = init_3d_array(nv,nv,nv);
  W = init_3d_array(nv,nv,nv);
  V = init_3d_array(nv,nv,nv);

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gkj = Gjk = Gk ^ Gj;
	Gji = Gij = Gi ^ Gj;
	Gik = Gki = Gi ^ Gk;

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      ij = T2.params->rowidx[I][J];
	      ji = T2.params->rowidx[J][I];
	      ik = T2.params->rowidx[I][K];
	      ki = T2.params->rowidx[K][I];
	      jk = T2.params->rowidx[J][K];
	      kj = T2.params->rowidx[K][J];

	      dijk = 0.0;
	      if(fIJ.params->rowtot[Gi])
		dijk += fIJ.matrix[Gi][i][i];
	      if(fIJ.params->rowtot[Gj])
		dijk += fIJ.matrix[Gj][j][j];
	      if(fIJ.params->rowtot[Gk])
		dijk += fIJ.matrix[Gk][k][k];

	      chksum = 0.0;
	      for(Ga=0; Ga < nirreps; Ga++) {
		for(Gb=0; Gb < nirreps; Gb++) {
		  Gc = Gi ^ Gj ^ Gk ^ Ga ^ Gb;

		  Gab = Gba = Ga ^ Gb;
		  Gac = Gca = Ga ^ Gc;
		  Gbc = Gcb = Gb ^ Gc;

		  for(a=0; a < virtpi[Ga]; a++) {
		    A = vir_off[Ga] + a;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			ab = Fints.params->colidx[A][B];
			ba = Fints.params->colidx[B][A];
			ac = Fints.params->colidx[A][C];
			ca = Fints.params->colidx[C][A];
			bc = Fints.params->colidx[B][C];
			cb = Fints.params->colidx[C][B];

			value = 0.0;

			/** Compute W and V **/

			/* +F_idab * t_kjcd */
			Gd = Gi ^ Ga ^ Gb;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  id = Fints.params->rowidx[I][D];
			  cd = T2.params->colidx[C][D];

			  F_idab = t_kjcd = 0.0;

			  if(Fints.params->rowtot[Gab] && Fints.params->coltot[Gab])
			    F_idab = Fints.matrix[Gab][id][ab];

			  if(T2.params->rowtot[Gkj] && T2.params->coltot[Gkj])
			    t_kjcd = T2.matrix[Gkj][kj][cd];

			  value += F_idab * t_kjcd;

			}

			/* +F_idac * t_jkbd */
			Gd = Gi ^ Ga ^ Gc;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  id = Fints.params->rowidx[I][D];
			  bd = T2.params->colidx[B][D];

			  F_idac = t_jkbd = 0.0;

			  if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_idac = Fints.matrix[Gac][id][ac];

			  if(T2.params->rowtot[Gjk] && T2.params->coltot[Gjk])
			    t_jkbd = T2.matrix[Gjk][jk][bd];

			  value += F_idac * t_jkbd;
			}

			/* +F_kdca * t_jibd */
			Gd = Gk ^ Gc ^ Ga;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  kd = Fints.params->rowidx[K][D];
			  bd = T2.params->colidx[B][D];

			  F_kdca = t_jibd = 0.0;

			  if(Fints.params->rowtot[Gca] && Fints.params->coltot[Gca])
			    F_kdca = Fints.matrix[Gca][kd][ca];

			  if(T2.params->rowtot[Gji] && T2.params->coltot[Gji])
			    t_jibd = T2.matrix[Gji][ji][bd];

			  value += F_kdca * t_jibd;
			}

			/* +F_kdcb * t_ijad */
			Gd = Gk ^ Gc ^ Gb;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  kd = Fints.params->rowidx[K][D];
			  ad = T2.params->colidx[A][D];

			  F_kdcb = t_ijad = 0.0;

			  if(Fints.params->rowtot[Gcb] && Fints.params->coltot[Gcb])
			    F_kdcb = Fints.matrix[Gcb][kd][cb];

			  if(T2.params->rowtot[Gij] && T2.params->coltot[Gij])
			    t_ijad = T2.matrix[Gij][ij][ad];

			  value += F_kdcb * t_ijad;
			}

			/* +F_jdbc * t_ikad */
			Gd = Gj ^ Gb ^ Gc;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  jd = Fints.params->rowidx[J][D];
			  ad = T2.params->colidx[A][D];

			  F_jdbc = t_ikad = 0.0;

			  if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_jdbc = Fints.matrix[Gbc][jd][bc];

			  if(T2.params->rowtot[Gik] && T2.params->coltot[Gik])
			    t_ikad = T2.matrix[Gik][ik][ad];

			  value += F_jdbc * t_ikad;
			}

			/* +F_jdba * t_kicd */
			Gd = Gj ^ Gb ^ Ga;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  jd = Fints.params->rowidx[J][D];
			  cd = T2.params->colidx[C][D];

			  F_jdba = t_kicd = 0.0;

			  if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_jdba = Fints.matrix[Gba][jd][ba];

			  if(T2.params->rowtot[Gki] && T2.params->coltot[Gki])
			    t_kicd = T2.matrix[Gki][ki][cd];

			  value += F_jdba * t_kicd;
			}

			/* -E_jklc * t_ilab */
			Gl = Gj ^ Gk ^ Gc;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  il = T2.params->rowidx[I][L];
			  lc = Eints.params->colidx[L][C];

			  E_jklc = t_ilab = 0.0;

			  if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_jklc = Eints.matrix[Gjk][jk][lc];

			  if(T2.params->rowtot[Gab] && T2.params->coltot[Gab])
			    t_ilab = T2.matrix[Gab][il][ab];

			  value -= E_jklc * t_ilab;
			}

			/* -E_kjlb * t_ilac */
			Gl = Gk ^ Gj ^ Gb;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  il = T2.params->rowidx[I][L];
			  lb = Eints.params->colidx[L][B];

			  E_kjlb = t_ilac = 0.0;

			  if(Eints.params->rowtot[Gkj] && Eints.params->coltot[Gkj])
			    E_kjlb = Eints.matrix[Gkj][kj][lb];

			  if(T2.params->rowtot[Gac] && T2.params->coltot[Gac])
			    t_ilac = T2.matrix[Gac][il][ac];

			  value -= E_kjlb * t_ilac;
			}

			/* -E_ijlb * t_klca */
			Gl = Gi ^ Gj ^ Gb;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  kl = T2.params->rowidx[K][L];
			  lb = Eints.params->colidx[L][B];

			  E_ijlb = t_klca = 0.0;

			  if(Eints.params->rowtot[Gij] && Eints.params->coltot[Gij])
			    E_ijlb = Eints.matrix[Gij][ij][lb];

			  if(T2.params->rowtot[Gca] && T2.params->coltot[Gca])
			    t_klca = T2.matrix[Gca][kl][ca];

			  value -= E_ijlb * t_klca;
			}

			/* -E_jila * t_klcb */
			Gl = Gj ^ Gi ^ Ga;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  kl = T2.params->rowidx[K][L];
			  la = Eints.params->colidx[L][A];

			  E_jila = t_klcb = 0.0;

			  if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_jila = Eints.matrix[Gji][ji][la];

			  if(T2.params->rowtot[Gcb] && T2.params->coltot[Gcb])
			    t_klcb = T2.matrix[Gcb][kl][cb];

			  value -= E_jila * t_klcb;
			}

			/* -E_kila * t_jlbc */
			Gl = Gk ^ Gi ^ Ga;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  jl = T2.params->rowidx[J][L];
			  la = Eints.params->colidx[L][A];

			  E_kila = t_jlbc = 0.0;

			  if(Eints.params->rowtot[Gki] && Eints.params->coltot[Gki])
			    E_kila = Eints.matrix[Gki][ki][la];

			  if(T2.params->rowtot[Gbc] && T2.params->coltot[Gbc])
			    t_jlbc = T2.matrix[Gbc][jl][bc];

			  value -= E_kila * t_jlbc;
			}

			/* -E_iklc * t_jlba */
			Gl = Gi ^ Gk ^ Gc;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  jl = T2.params->rowidx[J][L];
			  lc = Eints.params->colidx[L][C];

			  E_iklc = t_jlba = 0.0;

			  if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_iklc = Eints.matrix[Gik][ik][lc];

			  if(T2.params->rowtot[Gba] && T2.params->coltot[Gba])
			    t_jlba = T2.matrix[Gba][jl][ba];

			  value -= E_iklc * t_jlba;
			}

			W[A][B][C] = value;

			/* +t_ia * D_jkbc */
			if(Gi == Ga && Gjk == Gbc) {
			  t_ia = D_jkbc = 0.0;

			  if(T1.params->rowtot[Gi] && T1.params->coltot[Gi])
			    t_ia = T1.matrix[Gi][i][a];

			  if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    D_jkbc = Dints.matrix[Gjk][jk][bc];

			  value += t_ia * D_jkbc;

			}

			/* +t_jb * D_ikac */
			if(Gj == Gb && Gik == Gac) {
			  t_jb = D_ikac = 0.0;

			  if(T1.params->rowtot[Gj] && T1.params->coltot[Gj])
			    t_jb = T1.matrix[Gj][j][b];

			  if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikac = Dints.matrix[Gik][ik][ac];

			  value += t_jb * D_ikac;
			}

			/* +t_kc * D_ijab */
			if(Gk == Gc && Gij == Gab) {
			  t_kc = D_ijab = 0.0;

			  if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			    t_kc = T1.matrix[Gk][k][c];

			  if(Dints.params->rowtot[Gij] && Dints.params->coltot[Gij])
			    D_ijab = Dints.matrix[Gij][ij][ab];

			  value += t_kc * D_ijab;
			}

			V[A][B][C] = value;

		      } /* c */
		    } /* b */
		  } /* a */

		} /* Gb */
	      } /* Ga */

	      for(Ga=0; Ga < nirreps; Ga++) {
		for(Gb=0; Gb < nirreps; Gb++) {
		  Gc = Gi ^ Gj ^ Gk ^ Ga ^ Gb;

		  for(a=0; a < virtpi[Ga]; a++) {
		    A = vir_off[Ga] + a;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			/** Compute ET contributions **/

			value1 = (4.0 * W[A][B][C] + W[B][C][A] + W[C][A][B]);

			value2 = (V[A][B][C] - V[C][B][A]);

			denom = dijk;
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc])
			  denom -= fAB.matrix[Gc][c][c];

			ET += value1 * value2 / denom;

		      } /* c */
		    } /* b */
		  } /* a */

		} /* Gb */
	      } /* Ga */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  ET /= 3.0;

  free_3d_array(W, nv, nv);
  free_3d_array(V, nv, nv);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2, h);
    dpd_buf4_mat_irrep_close(&Fints, h);
    dpd_buf4_mat_irrep_close(&Eints, h);
    dpd_buf4_mat_irrep_close(&Dints, h);
  }

  dpd_buf4_close(&T2);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&Dints);

  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  return ET;
}
