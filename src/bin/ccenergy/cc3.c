#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

double ***init_3d_array(int p, int q, int r);
void free_3d_array(double ***A, int p, int q);

void cc3(void)
{
  int h, nirreps;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int Gi, Gj, Gk, Gl, Ga, Gb, Gc, Gd;
  int i, j, k, l, a, b, c, d;
  int I, J, K, L, A, B, C, D;
  int kj, jk, ji, ij, ik, ki;
  int Gkj, Gjk, Gji, Gij, Gik, Gki;
  int Gijk;
  int ab, ba, ac, ca, bc, cb;
  int Gab, Gba, Gac, Gca, Gbc, Gcb;
  int id, jd, kd, ad, bd, cd;
  int il, jl, kl, la, lb, lc, li, lk;
  int da, di, dj, dk;
  int Gad, Gdi, Gdj, Gdk, Glc, Gli, Glk;
  double value, F_val, t_val, E_val;
  double dijk, denom;
  double value_ia, value_ka, denom_ia, denom_ka;
  dpdfile2 fIJ, fAB, t1, t1new, Fme;
  dpdbuf4 T2, T2new, F, E, Dints, Wamef, Wmnie;
  double t2norm, t3norm;
  double ***T3;
  int nv;

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

  dpd_file2_init(&t1new, CC_OEI, 0, 0, 1, "CC3 tIA");  /* T3->T1 increment */
  dpd_file2_mat_init(&t1new);
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

  dpd_buf4_init(&T2new, CC_MISC, 0, 0, 5, 0, 5, 0, "CC3 tIjAb");
  dpd_buf4_scm(&T2new, 0);
  dpd_buf4_init(&F, CC_MISC, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (Ie,Ab)");
  dpd_buf4_init(&E, CC_MISC, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&Wamef, CC_MISC, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
  dpd_buf4_init(&Wmnie, CC_MISC, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2, h);
    dpd_buf4_mat_irrep_rd(&T2, h);

    dpd_buf4_mat_irrep_init(&F, h);
    dpd_buf4_mat_irrep_rd(&F, h);

    dpd_buf4_mat_irrep_init(&E, h);
    dpd_buf4_mat_irrep_rd(&E, h);

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);

    dpd_buf4_mat_irrep_init(&Wamef, h);
    dpd_buf4_mat_irrep_rd(&Wamef, h);

    dpd_buf4_mat_irrep_init(&Wmnie, h);
    dpd_buf4_mat_irrep_rd(&Wmnie, h);

    dpd_buf4_mat_irrep_init(&T2new, h);
  }

  for(h=0,nv=0; h < nirreps; h++) nv += virtpi[h];
  T3 = init_3d_array(nv, nv, nv);

  t3norm = 0.0;
  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gkj = Gjk = Gk ^ Gj;
	Gji = Gij = Gi ^ Gj;
	Gik = Gki = Gi ^ Gk;

	Gijk = Gi ^ Gj ^ Gk;

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

			denom = dijk;
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc])
			  denom -= fAB.matrix[Gc][c][c];

			ab = F.params->colidx[A][B];
			ba = F.params->colidx[B][A];
			ac = F.params->colidx[A][C];
			ca = F.params->colidx[C][A];
			bc = F.params->colidx[B][C];
			cb = F.params->colidx[C][B];

			value = 0.0;

			/* +F_idab * t_kjcd */
			Gd = Gi ^ Ga ^ Gb;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  id = F.params->rowidx[I][D];
			  cd = T2.params->colidx[C][D];

			  F_val = t_val = 0.0;

			  if(F.params->rowtot[Gab] && F.params->coltot[Gab])
			    F_val = F.matrix[Gab][id][ab];

			  if(T2.params->rowtot[Gkj] && T2.params->coltot[Gkj])
			    t_val = T2.matrix[Gkj][kj][cd];

			  value += F_val * t_val;

			}

			/* -E_jklc * t_ilab */
			Gl = Gj ^ Gk ^ Gc;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  il = T2.params->rowidx[I][L];
			  lc = E.params->colidx[L][C];

			  E_val = t_val = 0.0;

			  if(E.params->rowtot[Gjk] && E.params->coltot[Gjk])
			    E_val = E.matrix[Gjk][jk][lc];

			  if(T2.params->rowtot[Gab] && T2.params->coltot[Gab])
			    t_val = T2.matrix[Gab][il][ab];

			  value -= E_val * t_val;
			}

			/* +F_idac * t_jkbd */
			Gd = Gi ^ Ga ^ Gc;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  id = F.params->rowidx[I][D];
			  bd = T2.params->colidx[B][D];

			  F_val = t_val = 0.0;

			  if(F.params->rowtot[Gac] && F.params->coltot[Gac])
			    F_val = F.matrix[Gac][id][ac];

			  if(T2.params->rowtot[Gjk] && T2.params->coltot[Gjk])
			    t_val = T2.matrix[Gjk][jk][bd];

			  value += F_val * t_val;
			}

			/* -E_kjlb * t_ilac */
			Gl = Gk ^ Gj ^ Gb;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  il = T2.params->rowidx[I][L];
			  lb = E.params->colidx[L][B];

			  E_val = t_val = 0.0;

			  if(E.params->rowtot[Gkj] && E.params->coltot[Gkj])
			    E_val = E.matrix[Gkj][kj][lb];

			  if(T2.params->rowtot[Gac] && T2.params->coltot[Gac])
			    t_val = T2.matrix[Gac][il][ac];

			  value -= E_val * t_val;
			}

			/* +F_kdca * t_jibd */
			Gd = Gk ^ Gc ^ Ga;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  kd = F.params->rowidx[K][D];
			  bd = T2.params->colidx[B][D];

			  F_val = t_val = 0.0;

			  if(F.params->rowtot[Gca] && F.params->coltot[Gca])
			    F_val = F.matrix[Gca][kd][ca];

			  if(T2.params->rowtot[Gji] && T2.params->coltot[Gji])
			    t_val = T2.matrix[Gji][ji][bd];

			  value += F_val * t_val;
			}

			/* -E_ijlb * t_klca */
			Gl = Gi ^ Gj ^ Gb;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  kl = T2.params->rowidx[K][L];
			  lb = E.params->colidx[L][B];

			  E_val = t_val = 0.0;

			  if(E.params->rowtot[Gij] && E.params->coltot[Gij])
			    E_val = E.matrix[Gij][ij][lb];

			  if(T2.params->rowtot[Gca] && T2.params->coltot[Gca])
			    t_val = T2.matrix[Gca][kl][ca];

			  value -= E_val * t_val;
			}

			/* +F_kdcb * t_ijad */
			Gd = Gk ^ Gc ^ Gb;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  kd = F.params->rowidx[K][D];
			  ad = T2.params->colidx[A][D];

			  F_val = t_val = 0.0;

			  if(F.params->rowtot[Gcb] && F.params->coltot[Gcb])
			    F_val = F.matrix[Gcb][kd][cb];

			  if(T2.params->rowtot[Gij] && T2.params->coltot[Gij])
			    t_val = T2.matrix[Gij][ij][ad];

			  value += F_val * t_val;
			}

			/* -E_jila * t_klcb */
			Gl = Gj ^ Gi ^ Ga;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  kl = T2.params->rowidx[K][L];
			  la = E.params->colidx[L][A];

			  E_val = t_val = 0.0;

			  if(E.params->rowtot[Gji] && E.params->coltot[Gji])
			    E_val = E.matrix[Gji][ji][la];

			  if(T2.params->rowtot[Gcb] && T2.params->coltot[Gcb])
			    t_val = T2.matrix[Gcb][kl][cb];

			  value -= E_val * t_val;
			}

			/* +F_jdbc * t_ikad */
			Gd = Gj ^ Gb ^ Gc;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  jd = F.params->rowidx[J][D];
			  ad = T2.params->colidx[A][D];

			  F_val = t_val = 0.0;

			  if(F.params->rowtot[Gbc] && F.params->coltot[Gbc])
			    F_val = F.matrix[Gbc][jd][bc];

			  if(T2.params->rowtot[Gik] && T2.params->coltot[Gik])
			    t_val = T2.matrix[Gik][ik][ad];

			  value += F_val * t_val;
			}

			/* -E_kila * t_jlbc */
			Gl = Gk ^ Gi ^ Ga;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  jl = T2.params->rowidx[J][L];
			  la = E.params->colidx[L][A];

			  E_val = t_val = 0.0;

			  if(E.params->rowtot[Gki] && E.params->coltot[Gki])
			    E_val = E.matrix[Gki][ki][la];

			  if(T2.params->rowtot[Gbc] && T2.params->coltot[Gbc])
			    t_val = T2.matrix[Gbc][jl][bc];

			  value -= E_val * t_val;
			}

			/* +F_jdba * t_kicd */
			Gd = Gj ^ Gb ^ Ga;
			for(d=0; d < virtpi[Gd]; d++) {
			  D = vir_off[Gd] + d;

			  jd = F.params->rowidx[J][D];
			  cd = T2.params->colidx[C][D];

			  F_val = t_val = 0.0;

			  if(F.params->rowtot[Gba] && F.params->coltot[Gba])
			    F_val = F.matrix[Gba][jd][ba];

			  if(T2.params->rowtot[Gki] && T2.params->coltot[Gki])
			    t_val = T2.matrix[Gki][ki][cd];

			  value += F_val * t_val;
			}

			/* -E_iklc * t_jlba */
			Gl = Gi ^ Gk ^ Gc;
			for(l=0; l < occpi[Gl]; l++) {
			  L = occ_off[Gl] + l;

			  jl = T2.params->rowidx[J][L];
			  lc = E.params->colidx[L][C];

			  E_val = t_val = 0.0;

			  if(E.params->rowtot[Gik] && E.params->coltot[Gik])
			    E_val = E.matrix[Gik][ik][lc];

			  if(T2.params->rowtot[Gba] && T2.params->coltot[Gba])
			    t_val = T2.matrix[Gba][jl][ba];

			  value -= E_val * t_val;
			}

			value /= denom;
			t3norm += value * value;

			T3[A][B][C] = value;

		      } /* c */
		    } /* b */
		  } /* a */

		} /* Gb */
	      } /* Ga */

	      /*** T3 --> T1 Contributions ***/

	      for(Ga=0; Ga < nirreps; Ga++) {
		for(a=0; a < virtpi[Ga]; a++) {
		  A = vir_off[Ga] + a;

		  value_ia = 0.0;
		  value_ka = 0.0;
		  for(Gb=0; Gb < nirreps; Gb++) {
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;

		      Gc = Gijk ^ Ga ^ Gb;
		      Gbc = Gb ^ Gc;

		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			bc = Dints.params->colidx[B][C];

			if(Gi == Ga && Gjk == Gbc) {

			  if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    value_ia += T3[A][B][C] * Dints.matrix[Gjk][jk][bc];
			}

			if(Gk == Ga && Gji == Gbc) {

			  if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
			    value_ka -= T3[A][B][C] * Dints.matrix[Gji][ji][bc];
			}

		      } /* c */
		    } /* b */
		  } /* Gb */

		  denom_ia = denom_ka = 0.0;
		  if(fIJ.params->rowtot[Gi])
		    denom_ia += fIJ.matrix[Gi][i][i];

		  if(fIJ.params->rowtot[Gk])
		    denom_ka += fIJ.matrix[Gk][k][k];

		  if(fAB.params->rowtot[Ga]) {
		    denom_ia -= fAB.matrix[Ga][a][a];
		    denom_ka -= fAB.matrix[Ga][a][a];
		  }
		  value_ia /= denom_ia;
		  value_ka /= denom_ka;

		  if(t1new.params->rowtot[Gi] && t1new.params->coltot[Gi])
		    t1new.matrix[Gi][i][a] += value_ia;
		  if(t1new.params->rowtot[Gk] && t1new.params->coltot[Gk])
		    t1new.matrix[Gk][k][a] += value_ka;

		} /* a */
	      } /* Ga */

	      /*** T3 --> T2 Contributions ***/

	      for(Ga=0; Ga < nirreps; Ga++) {
	       	for(Gb=0; Gb < nirreps; Gb++) {
		  Gc = Gijk ^ Ga ^ Gb;
		  Gab = Ga ^ Gb;

		  for(a=0; a < virtpi[Ga]; a++) {
		    A = vir_off[Ga] + a;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;

                      ab = T2new.params->colidx[A][B];
                      ba = T2new.params->colidx[B][A];

		      if(Gij == Gab && Gk == Gc) { 

			denom = 0.0;
		       	if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
		       	if(fIJ.params->rowtot[Gj])
			  denom += fIJ.matrix[Gj][j][j];
		       	if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
		       	if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];

			value = 0.0;
		       	for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
			  if(Fme.params->rowtot[Gk] && Fme.params->coltot[Gc])
			    value += T3[A][B][C] * Fme.matrix[Gk][k][c];
			}
			value /= denom;

			if(T2new.params->rowtot[Gij] && T2new.params->coltot[Gab]) {
			  T2new.matrix[Gij][ij][ab] += value;
			  T2new.matrix[Gij][ji][ba] += value;
			}
		      }

                      if(Gjk == Gab && Gi == Gc) {

                        denom = 0.0;
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
                        if(fIJ.params->rowtot[Gj])
                          denom += fIJ.matrix[Gj][j][j];
                        if(fAB.params->rowtot[Ga])
                          denom -= fAB.matrix[Ga][a][a];
                        if(fAB.params->rowtot[Gb])
                          denom -= fAB.matrix[Gb][b][b];

                        value = 0.0;
                        for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
                          if(Fme.params->rowtot[Gi] && Fme.params->coltot[Gc])
                            value -= T3[A][B][C] * Fme.matrix[Gi][i][c];
			}
                        value /= denom;

                        if(T2new.params->rowtot[Gjk] && T2new.params->coltot[Gab]) {
                          T2new.matrix[Gjk][kj][ab] += value;
                          T2new.matrix[Gjk][jk][ba] += value;
                        }
                      }

		    } /* b */
		  } /* a */
	 	} /* Gb */
	      } /* Ga */

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gdi = Gd ^ Gi;
		Gdj = Gd ^ Gj;
		Gdk = Gd ^ Gk;
	       	for(Ga=0; Ga < nirreps; Ga++) {
		  Gad = Ga ^ Gd;

		  for(d=0; d < virtpi[Gd]; d++) {
		    D = vir_off[Gd] + d;

		    di = Wamef.params->rowidx[D][I];
		    dj = Wamef.params->rowidx[D][J];
		    dk = Wamef.params->rowidx[D][K];

		    for(a=0; a < virtpi[Ga]; a++) {
		      A = vir_off[Ga] + a;

		      ad = T2new.params->colidx[A][D];
		      da = T2new.params->colidx[D][A];

		      if(Gij == Gad) {
			value = 0.0;
			for(Gb=0; Gb < nirreps; Gb++) {
			  Gc = Gijk ^ Ga ^ Gb;
			  Gbc = Gb ^ Gc;
			  if(Gdk == Gbc) {
			    for(b=0; b < virtpi[Gb]; b++) {
			      B = vir_off[Gb] + b;
			      for(c=0; c < virtpi[Gc]; c++) {
			       	C = vir_off[Gc] + c;

				bc = Wamef.params->colidx[B][C];

				if(Wamef.params->rowtot[Gdk] && Wamef.params->coltot[Gdk])
				  value += 2.0 * T3[A][B][C] * Wamef.matrix[Gdk][dk][bc];
			      }
			    }
			  }
			}

			denom = 0.0;
			if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
			if(fIJ.params->rowtot[Gj])
			  denom += fIJ.matrix[Gj][j][j];
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gd])
			  denom -= fAB.matrix[Gd][d][d];
			value /= denom;

			if(T2new.params->rowtot[Gij] && T2new.params->coltot[Gij]) {
			  T2new.matrix[Gij][ij][ad] += value;
			  T2new.matrix[Gij][ji][da] += value;
			}
		      }

		      if(Gjk == Gad) {
                        value = 0.0; 
                        for(Gb=0; Gb < nirreps; Gb++) {
                          Gc = Gijk ^ Ga ^ Gb;
                          Gbc = Gb ^ Gc; 
                          if(Gdi == Gbc) {
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
                              for(c=0; c < virtpi[Gc]; c++) {
                                C = vir_off[Gc] + c;

				bc = Wamef.params->colidx[B][C];

                                if(Wamef.params->rowtot[Gdi] && Wamef.params->coltot[Gdi])
                                  value -= T3[A][B][C] * Wamef.matrix[Gdi][di][bc];
                              }   
                            } 
                          } 
                        } 
                        
                        denom = 0.0;
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
                        if(fIJ.params->rowtot[Gj])
                          denom += fIJ.matrix[Gj][j][j];
                        if(fAB.params->rowtot[Ga])
                          denom -= fAB.matrix[Ga][a][a];
                        if(fAB.params->rowtot[Gd])
                          denom -= fAB.matrix[Gd][d][d];
                        value /= denom;
                        
                        if(T2new.params->rowtot[Gjk] && T2new.params->coltot[Gjk]) {
                          T2new.matrix[Gjk][kj][ad] += value;
                          T2new.matrix[Gjk][jk][da] += value;
                        } 

		      }

		      if(Gik == Gad) {
                        value = 0.0;
                        for(Gb=0; Gb < nirreps; Gb++) {
                          Gc = Gijk ^ Ga ^ Gb;
                          Gbc = Gb ^ Gc;
                          if(Gdj == Gbc) {
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
                              for(c=0; c < virtpi[Gc]; c++) {
                                C = vir_off[Gc] + c;

				bc = Wamef.params->colidx[B][C];

                                if(Wamef.params->rowtot[Gdj] && Wamef.params->coltot[Gdj])
                                  value -= T3[A][B][C] * Wamef.matrix[Gdj][dj][bc];
                              }
                            }
                          }
                        }

                        denom = 0.0;
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
                        if(fIJ.params->rowtot[Gi])
                          denom += fIJ.matrix[Gi][i][i];
                        if(fAB.params->rowtot[Ga])
                          denom -= fAB.matrix[Ga][a][a];
                        if(fAB.params->rowtot[Gd])
                          denom -= fAB.matrix[Gd][d][d];
                        value /= denom;

                        if(T2new.params->rowtot[Gik] && T2new.params->coltot[Gik]) {
                          T2new.matrix[Gki][ki][da] += value;
                          T2new.matrix[Gki][ik][ad] += value;
                        }
		      }

		    } /* a */
		  } /* d */
		} /* Ga */
	      } /* Gd */

	      for(Gl=0; Gl < nirreps; Gl++) {
		Gli = Gl ^ Gi;
		Glk = Gl ^ Gk;

		for(Ga=0; Ga < nirreps; Ga++) {
		  for(Gb=0; Gb < nirreps; Gb++) {
		    Gab = Ga ^ Gb;

		    if(Gli == Gab) {

		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;

			li = T2new.params->rowidx[L][I];
			il = T2new.params->rowidx[I][L];

			for(a=0; a < virtpi[Ga]; a++) {
			  A = vir_off[Ga] + a;
			  for(b=0; b < virtpi[Gb]; b++) {
			    B = vir_off[Gb] + b;

			    ab = T2new.params->colidx[A][B];
			    ba = T2new.params->colidx[B][A];

			    denom = 0.0;
			    if(fIJ.params->rowtot[Gl])
			      denom += fIJ.matrix[Gl][l][l];
			    if(fIJ.params->rowtot[Gi])
			      denom += fIJ.matrix[Gi][i][i];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];

			    value = 0.0;
			    for(Gc=0; Gc < nirreps; Gc++) {
			      Glc = Gl ^ Gc;

			      if(Gjk == Glc) {

				for(c=0; c < virtpi[Gc]; c++) {
				  C = vir_off[Gc] + c;
				  lc = Wmnie.params->colidx[L][C];

				  if(Wmnie.params->rowtot[Gjk] && Wmnie.params->coltot[Gjk])
				    value -= T3[A][B][C] * Wmnie.matrix[Gjk][jk][lc];
				}
			      }

			    } /* Gc */

			    value *= 2.0;
			    value /= denom;

			    if(T2new.params->rowtot[Gli] && T2new.params->coltot[Gli]) {
			      T2new.matrix[Gli][li][ba] += value;
			      T2new.matrix[Gli][il][ab] += value;
			    }

			  } /* b */
			} /* a */
		      } /* l */

		    } /* if(Gli == Gab) */

		    if(Glk == Gab) {

		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;

			lk = T2new.params->rowidx[L][K];
			kl = T2new.params->rowidx[K][L];

			for(a=0; a < virtpi[Ga]; a++) {
			  A = vir_off[Ga] + a;
			  for(b=0; b < virtpi[Gb]; b++) {
			    B = vir_off[Gb] + b;

			    ab = T2new.params->colidx[A][B];
			    ba = T2new.params->colidx[B][A];

			    denom = 0.0;
			    if(fIJ.params->rowtot[Gl])
			      denom += fIJ.matrix[Gl][l][l];
			    if(fIJ.params->rowtot[Gk])
			      denom += fIJ.matrix[Gk][k][k];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];

			    value = 0.0;
			    for(Gc=0; Gc < nirreps; Gc++) {
			      Glc = Gl ^ Gc;

			      if(Gji == Glc) {

				for(c=0; c < virtpi[Gc]; c++) {
				  C = vir_off[Gc] + c;
				  lc = Wmnie.params->colidx[L][C];

				  if(Wmnie.params->rowtot[Gji] && Wmnie.params->coltot[Gji])
				    value += T3[A][B][C] * Wmnie.matrix[Gji][ji][lc];
				}
			      }

			    } /* Gc */

			    value /= denom;

			    if(T2new.params->rowtot[Glk] && T2new.params->coltot[Glk]) {
			      T2new.matrix[Glk][lk][ba] += value;
			      T2new.matrix[Glk][kl][ab] += value;
			    }

			  } /* b */
			} /* a */
		      } /* l */

		    } /* if(Glk==Gab) */

		    if(Gli == Gab) {

		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;

			li = T2new.params->rowidx[L][I];
			il = T2new.params->rowidx[I][L];

			for(a=0; a < virtpi[Ga]; a++) {
			  A = vir_off[Ga] + a;
			  for(b=0; b < virtpi[Gb]; b++) {
			    B = vir_off[Gb] + b;

			    ab = T2new.params->colidx[A][B];
			    ba = T2new.params->colidx[B][A];

			    denom = 0.0;
			    if(fIJ.params->rowtot[Gl])
			      denom += fIJ.matrix[Gl][l][l];
			    if(fIJ.params->rowtot[Gi])
			      denom += fIJ.matrix[Gi][i][i];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];

			    value = 0.0;
			    for(Gc=0; Gc < nirreps; Gc++) {
			      Glc = Gl ^ Gc;

			      if(Gjk == Glc) {

				for(c=0; c < virtpi[Gc]; c++) {
				  C = vir_off[Gc] + c;
				  lc = Wmnie.params->colidx[L][C];

				  if(Wmnie.params->rowtot[Gjk] && Wmnie.params->coltot[Gjk])
				    value += T3[A][B][C] * Wmnie.matrix[Gjk][kj][lc];
				}
			      }

			    } /* Gc */

			    value /= denom;

			    if(T2new.params->rowtot[Gli] && T2new.params->coltot[Gli]) {
			      T2new.matrix[Gli][li][ba] += value;
			      T2new.matrix[Gli][il][ab] += value;
			    }

			  } /* b */
			} /* a */
		      } /* l */

		    } /* if(Gli==Gab) */

		  } /* Gb */
		} /* Ga */
	      } /* Gl */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  /*
  fprintf(outfile, "t3norm = %20.10f\n", sqrt(t3norm));
  */
  free_3d_array(T3, nv, nv);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2, h);
    dpd_buf4_mat_irrep_close(&E, h);
    dpd_buf4_mat_irrep_close(&F, h);
    dpd_buf4_mat_irrep_close(&Dints, h);
    dpd_buf4_mat_irrep_close(&Wamef, h);
    dpd_buf4_mat_irrep_close(&Wmnie, h);

    dpd_buf4_mat_irrep_wrt(&T2new, h);
    dpd_buf4_mat_irrep_close(&T2new, h);
  }

  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&E);
  dpd_buf4_close(&Dints);
  dpd_buf4_close(&Wamef);
  dpd_buf4_close(&Wmnie);

  /*
  fprintf(outfile, "t2norm = %20.10f\n", sqrt(dpd_buf4_dot_self(&T2new)));
  */

  /* Update the amplitudes */
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);

  dpd_file2_mat_wrt(&t1new);
  dpd_file2_mat_close(&t1new);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1new, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1new);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
}
