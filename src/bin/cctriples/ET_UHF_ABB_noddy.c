#include <stdio.h>
#include <math.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double ET_UHF_ABB_noddy(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Ge, Gm;
  int Gji, Gij, Gjk, Gik, Gbc, Gac, Gba, Gab;
  int I, J, K, A, B, C, E, M;
  int i, j, k, a, b, c, e, m;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int ae, be, eb, ce, ec, ie, je, ke;
  int im, jm, mj, km, mk, ma, mb, mc;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off;
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
  double value_c, value_d, denom, ET_ABB;
  double t_ijae, t_ijeb, t_ijec, t_jkbe, t_jkce, t_ikae, t_ikeb, t_ikec;
  double F_kebc, F_keca, F_keba, F_ieac, F_ieab, F_jebc, F_jeca, F_jeba;
  double t_imac, t_imab, t_jmbc, t_mjac, t_mjab, t_kmbc, t_mkac, t_mkab;
  double E_jkmb, E_jkmc, E_kima, E_ikmb, E_ikmc, E_jima, E_ijmb, E_ijmc;
  double t_ia, t_jb, t_jc, t_kb, t_kc;
  double D_jkbc, D_ikac, D_ikab, D_ijac, D_ijab;
  dpdbuf4 T2AB, T2BB;
  dpdbuf4 FBBints, FABints, FBAints;
  dpdbuf4 EBBints, EABints, EBAints;
  dpdbuf4 DBBints, DABints;
  dpdfile2 T1A, T1B, fIJ, fij, fAB, fab;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; 
  avirtpi = moinfo.avirtpi;
  aocc_off = moinfo.aocc_off;
  avir_off = moinfo.avir_off;
  boccpi = moinfo.boccpi; 
  bvirtpi = moinfo.bvirtpi;
  bocc_off = moinfo.bocc_off;
  bvir_off = moinfo.bvir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fij);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_rd(&fab);

  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1A);
  dpd_file2_mat_rd(&T1A);
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_mat_init(&T1B);
  dpd_file2_mat_rd(&T1B);

  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");

  dpd_buf4_init(&FBBints, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  dpd_buf4_init(&FABints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

  dpd_buf4_init(&EBBints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&EABints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

  dpd_buf4_init(&DBBints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  dpd_buf4_init(&DABints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BB, h);
    dpd_buf4_mat_irrep_rd(&T2BB, h);

    dpd_buf4_mat_irrep_init(&T2AB, h);
    dpd_buf4_mat_irrep_rd(&T2AB, h);

    dpd_buf4_mat_irrep_init(&FBBints, h);
    dpd_buf4_mat_irrep_rd(&FBBints, h);

    dpd_buf4_mat_irrep_init(&FABints, h);
    dpd_buf4_mat_irrep_rd(&FABints, h);

    dpd_buf4_mat_irrep_init(&FBAints, h);
    dpd_buf4_mat_irrep_rd(&FBAints, h);

    dpd_buf4_mat_irrep_init(&EBBints, h);
    dpd_buf4_mat_irrep_rd(&EBBints, h);

    dpd_buf4_mat_irrep_init(&EABints, h);
    dpd_buf4_mat_irrep_rd(&EABints, h);

    dpd_buf4_mat_irrep_init(&EBAints, h);
    dpd_buf4_mat_irrep_rd(&EBAints, h);

    dpd_buf4_mat_irrep_init(&DBBints, h);
    dpd_buf4_mat_irrep_rd(&DBBints, h);

    dpd_buf4_mat_irrep_init(&DABints, h);
    dpd_buf4_mat_irrep_rd(&DABints, h);
  }

  cnt = 0;
  ET_ABB = 0.0;

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {
	Gij = Gji = Gi ^ Gj;
	Gjk = Gj ^ Gk;
	Gik = Gi ^ Gk;

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < boccpi[Gj]; j++) {
	    J = bocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      ij = EABints.params->rowidx[I][J];
	      ji = EBAints.params->rowidx[J][I];
	      jk = EBBints.params->rowidx[J][K];
	      kj = EBBints.params->rowidx[K][J];
	      ik = EABints.params->rowidx[I][K];
	      ki = EBAints.params->rowidx[K][I];

	      for(Ga=0; Ga < nirreps; Ga++) {
		for(Gb=0; Gb < nirreps; Gb++) {
		  Gc = Gi ^ Gj ^ Gk ^ Ga ^ Gb;

		  Gbc = Gb^Gc;
		  Gac = Ga^Gc;
		  Gba = Gb^Ga;

		  for(a=0; a < avirtpi[Ga]; a++) {
		    A = avir_off[Ga] + a;
		    for(b=0; b < bvirtpi[Gb]; b++) {
		      B = bvir_off[Gb] + b;
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;

			ab = FABints.params->colidx[A][B];
			ba = FBAints.params->colidx[B][A];
			bc = FBBints.params->colidx[B][C];
			cb = FBBints.params->colidx[C][B];
			ac = FABints.params->colidx[A][C];
			ca = FBAints.params->colidx[C][A];

			value_c = 0.0;

			/** <ov||vv> --> connected triples **/

                        /* +t_jkbe * F_IeAc */
                        Ge = Gj ^ Gk ^ Gb;
                        for(e=0; e < bvirtpi[Ge]; e++) {
                          E = bvir_off[Ge] + e;

			  be = T2BB.params->colidx[B][E];
			  ie = FABints.params->rowidx[I][E];

                          t_jkbe = F_ieac = 0.0;

                          if(T2BB.params->rowtot[Gjk] && T2BB.params->coltot[Gjk])
                            t_jkbe = T2BB.matrix[Gjk][jk][be];

                          if(FABints.params->rowtot[Gac] && FABints.params->coltot[Gac])
			    F_ieac = FABints.matrix[Gac][ie][ac];
 
                          value_c += t_jkbe * F_ieac;
                        }

                        /* -t_jkce * F_IeAb */
                        Ge = Gj ^ Gk ^ Gc;
                        for(e=0; e < bvirtpi[Ge]; e++) {
                          E = bvir_off[Ge] + e;

			  ce = T2BB.params->colidx[C][E];
			  ie = FABints.params->rowidx[I][E];

                          t_jkce = F_ieab = 0.0;

                          if(T2BB.params->rowtot[Gjk] && T2BB.params->coltot[Gjk])
                            t_jkce = T2BB.matrix[Gjk][jk][ce];

                          if(FABints.params->rowtot[Gba] && FABints.params->coltot[Gba])
			    F_ieab = FABints.matrix[Gba][ie][ab];
 
                          value_c -= t_jkce * F_ieab;
                        }

                        /* +t_IkAe * F_jebc */
                        Ge = Gi ^ Gk ^ Ga;
                        for(e=0; e < bvirtpi[Ge]; e++) {
                          E = bvir_off[Ge] + e;

			  ae = T2AB.params->colidx[A][E];
			  je = FBBints.params->rowidx[J][E];

                          t_ikae = F_jebc = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikae = T2AB.matrix[Gik][ik][ae];

                          if(FBBints.params->rowtot[Gbc] && FBBints.params->coltot[Gbc])
			    F_jebc = FBBints.matrix[Gbc][je][bc];
 
                          value_c += t_ikae * F_jebc;
                        }

                        /* -t_IkEb * F_jEcA */
                        Ge = Gi ^ Gk ^ Gb;
                        for(e=0; e < avirtpi[Ge]; e++) {
                          E = avir_off[Ge] + e;

			  eb = T2AB.params->colidx[E][B];
			  je = FBAints.params->rowidx[J][E];

                          t_ikeb = F_jeca = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikeb = T2AB.matrix[Gik][ik][eb];

                          if(FBAints.params->rowtot[Gac] && FBAints.params->coltot[Gac])
			    F_jeca = FBAints.matrix[Gac][je][ca];
 
                          value_c -= t_ikeb * F_jeca;
                        }

                        /* +t_IkEc * F_jEbA */
                        Ge = Gi ^ Gk ^ Gc;
                        for(e=0; e < avirtpi[Ge]; e++) {
                          E = avir_off[Ge] + e;

			  ec = T2AB.params->colidx[E][C];
			  je = FBAints.params->rowidx[J][E];

                          t_ikec = F_jeba = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikec = T2AB.matrix[Gik][ik][ec];

                          if(FBAints.params->rowtot[Gba] && FBAints.params->coltot[Gba])
			    F_jeba = FBAints.matrix[Gba][je][ba];
 
                          value_c += t_ikec * F_jeba;
                        }

                        /* -t_IjAe * F_kebc */
                        Ge = Gi ^ Gj ^ Ga;
                        for(e=0; e < bvirtpi[Ge]; e++) {
                          E = bvir_off[Ge] + e;

			  ae = T2AB.params->colidx[A][E];
			  ke = FBBints.params->rowidx[K][E];

                          t_ijae = F_kebc = 0.0;

                          if(T2AB.params->rowtot[Gij] && T2AB.params->coltot[Gij])
                            t_ijae = T2AB.matrix[Gij][ij][ae];

                          if(FBBints.params->rowtot[Gbc] && FBBints.params->coltot[Gbc])
			    F_kebc = FBBints.matrix[Gbc][ke][bc];
 
                          value_c -= t_ijae * F_kebc;
                        }

                        /* +t_IjEb * F_kEcA */
                        Ge = Gi ^ Gj ^ Gb;
                        for(e=0; e < avirtpi[Ge]; e++) {
                          E = avir_off[Ge] + e;

			  eb = T2AB.params->colidx[E][B];
			  ke = FBAints.params->rowidx[K][E];

                          t_ijeb = F_keca = 0.0;

                          if(T2AB.params->rowtot[Gij] && T2AB.params->coltot[Gij])
                            t_ijeb = T2AB.matrix[Gij][ij][eb];

                          if(FBAints.params->rowtot[Gac] && FBAints.params->coltot[Gac])
			    F_keca = FBAints.matrix[Gac][ke][ca];
 
                          value_c += t_ijeb * F_keca;
                        }

                        /* -t_IjEc * F_kEbA */
                        Ge = Gi ^ Gj ^ Gc;
                        for(e=0; e < avirtpi[Ge]; e++) {
                          E = avir_off[Ge] + e;

			  ec = T2AB.params->colidx[E][C];
			  ke = FBAints.params->rowidx[K][E];

                          t_ijec = F_keba = 0.0;

                          if(T2AB.params->rowtot[Gij] && T2AB.params->coltot[Gij])
                            t_ijec = T2AB.matrix[Gij][ij][ec];

                          if(FBAints.params->rowtot[Gba] && FBAints.params->coltot[Gba])
			    F_keba = FBAints.matrix[Gba][ke][ba];
 
                          value_c -= t_ijec * F_keba;
                        }

			/** <oo||ov> --> connected triples **/

                        /* +t_ImAc * E_jkmb */
			Gm = Gi ^ Ga ^ Gc;
			for(m=0; m < boccpi[Gm]; m++) {
                          M = bocc_off[Gm] + m;

                          im = T2AB.params->rowidx[I][M];
                          mb = EBBints.params->colidx[M][B];

                          t_imac = E_jkmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_imac = T2AB.matrix[Gac][im][ac];

                          if(EBBints.params->rowtot[Gjk] && EBBints.params->coltot[Gjk])
			    E_jkmb = EBBints.matrix[Gjk][jk][mb];

			  value_c += t_imac * E_jkmb;
			}

                        /* -t_ImAb * E_jkmc */
			Gm = Gi ^ Gb ^ Ga;
			for(m=0; m < boccpi[Gm]; m++) {
                          M = bocc_off[Gm] + m;

                          im = T2AB.params->rowidx[I][M];
                          mc = EBBints.params->colidx[M][C];

                          t_imab = E_jkmc = 0.0;

                          if(T2AB.params->rowtot[Gba] && T2AB.params->coltot[Gba])
			    t_imab = T2AB.matrix[Gba][im][ab];

                          if(EBBints.params->rowtot[Gjk] && EBBints.params->coltot[Gjk])
			    E_jkmc = EBBints.matrix[Gjk][jk][mc];

                          value_c -= t_imab * E_jkmc;
			}

                        /* -t_jmbc * E_kImA */
			Gm = Gj ^ Gb ^ Gc;
			for(m=0; m < boccpi[Gm]; m++) {
                          M = bocc_off[Gm] + m;

                          jm = T2BB.params->rowidx[J][M];
                          ma = EBAints.params->colidx[M][A];

                          t_jmbc = E_kima = 0.0;

                          if(T2BB.params->rowtot[Gbc] && T2BB.params->coltot[Gbc])
			    t_jmbc = T2BB.matrix[Gbc][jm][bc];

                          if(EBAints.params->rowtot[Gik] && EBAints.params->coltot[Gik])
			    E_kima = EBAints.matrix[Gik][ki][ma];

                          value_c -= t_jmbc * E_kima;
			}

                        /* +t_MjAc * E_IkMb */
			Gm = Gj ^ Ga ^ Gc;
			for(m=0; m < aoccpi[Gm]; m++) {
                          M = aocc_off[Gm] + m;

                          mj = T2AB.params->rowidx[M][J];
                          mb = EABints.params->colidx[M][B];

                          t_mjac = E_ikmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_mjac = T2AB.matrix[Gac][mj][ac];

                          if(EABints.params->rowtot[Gik] && EABints.params->coltot[Gik])
			    E_ikmb = EABints.matrix[Gik][ik][mb];

                          value_c += t_mjac * E_ikmb;
			}

                        /* -t_MjAb * E_IkMc */
			Gm = Gj ^ Gb ^ Ga;
			for(m=0; m < aoccpi[Gm]; m++) {
                          M = aocc_off[Gm] + m;

                          mj = T2AB.params->rowidx[M][J];
                          mc = EABints.params->colidx[M][C];

                          t_mjab = E_ikmc = 0.0;

                          if(T2AB.params->rowtot[Gba] && T2AB.params->coltot[Gba])
			    t_mjab = T2AB.matrix[Gba][mj][ab];

                          if(EABints.params->rowtot[Gik] && EABints.params->coltot[Gik])
			    E_ikmc = EABints.matrix[Gik][ik][mc];

                          value_c -= t_mjab * E_ikmc;
			}

                        /* +t_kmbc * E_jImA */
			Gm = Gk ^ Gb ^ Gc;
			for(m=0; m < boccpi[Gm]; m++) {
                          M = bocc_off[Gm] + m;

                          km = T2BB.params->rowidx[K][M];
                          ma = EBAints.params->colidx[M][A];

                          t_kmbc = E_jima = 0.0;

                          if(T2BB.params->rowtot[Gbc] && T2BB.params->coltot[Gbc])
			    t_kmbc = T2BB.matrix[Gbc][km][bc];

                          if(EBAints.params->rowtot[Gji] && EBAints.params->coltot[Gji])
			    E_jima = EBAints.matrix[Gji][ji][ma];

                          value_c += t_kmbc * E_jima;
			}

                        /* -t_MkAc * E_IjMb */
			Gm = Gk ^ Ga ^ Gc;
			for(m=0; m < aoccpi[Gm]; m++) {
                          M = aocc_off[Gm] + m;

                          mk = T2AB.params->rowidx[M][K];
                          mb = EABints.params->colidx[M][B];

                          t_mkac = E_ijmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_mkac = T2AB.matrix[Gac][mk][ac];

                          if(EABints.params->rowtot[Gji] && EABints.params->coltot[Gji])
			    E_ijmb = EABints.matrix[Gji][ij][mb];

                          value_c -= t_mkac * E_ijmb;
			}

                        /* +t_MkAb * E_IjMc */
			Gm = Gk ^ Gb ^ Ga;
			for(m=0; m < aoccpi[Gm]; m++) {
                          M = aocc_off[Gm] + m;

                          mk = T2AB.params->rowidx[M][K];
                          mc = EABints.params->colidx[M][C];

                          t_mkab = E_ijmc = 0.0;

                          if(T2AB.params->rowtot[Gba] && T2AB.params->coltot[Gba])
			    t_mkab = T2AB.matrix[Gba][mk][ab];

                          if(EABints.params->rowtot[Gji] && EABints.params->coltot[Gji])
			    E_ijmc = EABints.matrix[Gji][ij][mc];

                          value_c += t_mkab * E_ijmc;
			}

			/** disconnected triples **/

			value_d = 0.0;

                        /* +t_IA * D_jkbc */
			if(Gi == Ga && Gjk == Gbc) {
			  t_ia = D_jkbc = 0.0;

                          if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			    t_ia = T1A.matrix[Gi][i][a];

                          if(DBBints.params->rowtot[Gjk] && DBBints.params->coltot[Gjk])
			    D_jkbc = DBBints.matrix[Gjk][jk][bc];

                          value_d += t_ia * D_jkbc;
			}

			/* +t_jb * D_IkAc */
			if(Gj == Gb && Gik == Gac) {
			  t_jb = D_ikac = 0.0;

                          if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj])
			    t_jb = T1B.matrix[Gj][j][b];

                          if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
			    D_ikac = DABints.matrix[Gik][ik][ac];

                          value_d += t_jb * D_ikac;
			}

			/* -t_jc * D_IkAb */
			if(Gj == Gc && Gik == Gba) {
			  t_jc = D_ikab = 0.0;

                          if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj])
			    t_jc = T1B.matrix[Gj][j][c];

                          if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
			    D_ikab = DABints.matrix[Gik][ik][ab];

                          value_d -= t_jc * D_ikab;
			}

			/* -t_kb * D_IjAc */
			if(Gk == Gb && Gji == Gac) {
			  t_kb = D_ijac = 0.0;

                          if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			    t_kb = T1B.matrix[Gk][k][b];

                          if(DABints.params->rowtot[Gji] && DABints.params->coltot[Gji])
			    D_ijac = DABints.matrix[Gji][ij][ac];

                          value_d -= t_kb * D_ijac;
			}

			/* +t_kc * D_IjAb */
			if(Gk == Gc && Gji == Gba) {
			  t_kc = D_ijab = 0.0;

                          if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			    t_kc = T1B.matrix[Gk][k][c];

                          if(DABints.params->rowtot[Gji] && DABints.params->coltot[Gji])
			    D_ijab = DABints.matrix[Gji][ij][ab];

                          value_d += t_kc * D_ijab;
			}

			/*
			  if(fabs(value_c) > 1e-7) {
			  cnt++;
			  fprintf(outfile, "%d %d %d %d %d %d %20.14f\n", I, J, K, A, B, C, value_c);
			  }
			*/

			/* Compute the Fock denominator */
			denom = 0.0;
			if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
			if(fij.params->rowtot[Gj])
			  denom += fij.matrix[Gj][j][j];
			if(fij.params->rowtot[Gk])
			  denom += fij.matrix[Gk][k][k];
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fab.params->rowtot[Gb])
			  denom -= fab.matrix[Gb][b][b];
			if(fab.params->rowtot[Gc])
			  denom -= fab.matrix[Gc][c][c];

			/*			ET_ABB += (value_d + value_c) * value_c / denom; */
			ET_ABB += (value_c + value_d) * value_c / denom;

		      } /* c */
		    } /* ab */
		  } /* Gab */

		} /* k */
	      } /* j */
	    } /* i */

	  } /* Gb */
	} /* Ga */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  /*  fprintf(outfile, "cnt = %d\n", cnt); */
  ET_ABB /= 4.0;

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2BB, h);
    dpd_buf4_mat_irrep_close(&T2AB, h);
    dpd_buf4_mat_irrep_close(&FBBints, h);
    dpd_buf4_mat_irrep_close(&FABints, h);
    dpd_buf4_mat_irrep_close(&FBAints, h);
    dpd_buf4_mat_irrep_close(&EBBints, h);
    dpd_buf4_mat_irrep_close(&EABints, h);
    dpd_buf4_mat_irrep_close(&EBAints, h);
    dpd_buf4_mat_irrep_close(&DBBints, h);
    dpd_buf4_mat_irrep_close(&DABints, h);
  }

  dpd_buf4_close(&T2BB);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&FBBints);
  dpd_buf4_close(&FABints);
  dpd_buf4_close(&FBAints);
  dpd_buf4_close(&EBBints);
  dpd_buf4_close(&EABints);
  dpd_buf4_close(&EBAints);
  dpd_buf4_close(&DBBints);
  dpd_buf4_close(&DABints);

  dpd_file2_mat_close(&T1A);
  dpd_file2_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1B);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);

  return ET_ABB;
}
