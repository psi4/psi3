/*! \file ET_UHF_AAB_noddy.c
    \ingroup (CCTRIPLES)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double ET_UHF_AAB_noddy(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Ge, Gm;
  int Gji, Gij, Gjk, Gik, Gbc, Gac, Gba;
  int I, J, K, A, B, C, E, M;
  int i, j, k, a, b, c, e, m;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int ae, be, ec, ke, ie, je;
  int im, jm, km, mk, ma, mb, mc;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off;
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
  double value_c, value_d, denom, ET_AAB;
  double t_ijae, t_ijbe, t_jkae, t_jkbe, t_jkec, t_ikae, t_ikbe, t_ikec;
  double F_kecb, F_keca, F_iebc, F_ieac, F_ieab, F_jebc, F_jeac, F_jeab;
  double t_imbc, t_imac, t_imba, t_jmbc, t_jmac, t_jmba, t_mkbc, t_mkac;
  double E_kjma, E_kjmb, E_jkmc, E_kima, E_kimb, E_ikmc, E_jima, E_jimb;
  double t_ia, t_ib, t_ja, t_jb, t_kc;
  double D_jkbc, D_jkac, D_ikbc, D_ikac, D_jiba;
  dpdbuf4 T2AB, T2AA;
  dpdbuf4 FAAints, FABints, FBAints;
  dpdbuf4 EAAints, EABints, EBAints;
  dpdbuf4 DAAints, DABints;
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

  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");

  dpd_buf4_init(&FAAints, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_buf4_init(&FABints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

  dpd_buf4_init(&EAAints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
  dpd_buf4_init(&EABints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

  dpd_buf4_init(&DAAints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  dpd_buf4_init(&DABints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2AA, h);
    dpd_buf4_mat_irrep_rd(&T2AA, h);

    dpd_buf4_mat_irrep_init(&T2AB, h);
    dpd_buf4_mat_irrep_rd(&T2AB, h);

    dpd_buf4_mat_irrep_init(&FAAints, h);
    dpd_buf4_mat_irrep_rd(&FAAints, h);

    dpd_buf4_mat_irrep_init(&FABints, h);
    dpd_buf4_mat_irrep_rd(&FABints, h);

    dpd_buf4_mat_irrep_init(&FBAints, h);
    dpd_buf4_mat_irrep_rd(&FBAints, h);

    dpd_buf4_mat_irrep_init(&EAAints, h);
    dpd_buf4_mat_irrep_rd(&EAAints, h);

    dpd_buf4_mat_irrep_init(&EABints, h);
    dpd_buf4_mat_irrep_rd(&EABints, h);

    dpd_buf4_mat_irrep_init(&EBAints, h);
    dpd_buf4_mat_irrep_rd(&EBAints, h);

    dpd_buf4_mat_irrep_init(&DAAints, h);
    dpd_buf4_mat_irrep_rd(&DAAints, h);

    dpd_buf4_mat_irrep_init(&DABints, h);
    dpd_buf4_mat_irrep_rd(&DABints, h);
  }

  cnt = 0;
  ET_AAB = 0.0;

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gij = Gji = Gi ^ Gj;
	Gjk = Gj ^ Gk;
	Gik = Gi ^ Gk;

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < aoccpi[Gj]; j++) {
	    J = aocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      if(I >= J) {

		ij = EAAints.params->rowidx[I][J];
		ji = EAAints.params->rowidx[J][I];
		jk = EABints.params->rowidx[J][K];
		kj = EBAints.params->rowidx[K][J];
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
		      for(b=0; b < avirtpi[Gb]; b++) {
			B = avir_off[Gb] + b;
			for(c=0; c < bvirtpi[Gc]; c++) {
			  C = bvir_off[Gc] + c;

			  /*			  if(A >= B) { */

			    ab = FAAints.params->colidx[A][B];
			    ba = FAAints.params->colidx[B][A];
			    bc = FABints.params->colidx[B][C];
			    cb = FBAints.params->colidx[C][B];
			    ac = FABints.params->colidx[A][C];
			    ca = FBAints.params->colidx[C][A];

			    value_c = 0.0;

			    /** <ov||vv> --> connected triples **/

			    /* -t_JkAe * F_IeBc */
			    Ge = Gj ^ Gk ^ Ga;
			    for(e=0; e < bvirtpi[Ge]; e++) {
			      E = bvir_off[Ge] + e;

			      ae = T2AB.params->colidx[A][E];
			      ie = FABints.params->rowidx[I][E];

			      t_jkae = F_iebc = 0.0;

			      if(T2AB.params->rowtot[Gjk] && T2AB.params->coltot[Gjk])
				t_jkae = T2AB.matrix[Gjk][jk][ae];

			      if(FABints.params->rowtot[Gbc] && FABints.params->coltot[Gbc])
				F_iebc = FABints.matrix[Gbc][ie][bc];
 
			      value_c -= t_jkae * F_iebc;
			    }

			    /* +t_JkBe * F_IeAc */
			    Ge = Gj ^ Gk ^ Gb;
			    for(e=0; e < bvirtpi[Ge]; e++) {
			      E = bvir_off[Ge] + e;

			      be = T2AB.params->colidx[B][E];
			      ie = FABints.params->rowidx[I][E];

			      t_jkbe = F_ieac = 0.0;

			      if(T2AB.params->rowtot[Gjk] && T2AB.params->coltot[Gjk])
				t_jkbe = T2AB.matrix[Gjk][jk][be];

			      if(FABints.params->rowtot[Gac] && FABints.params->coltot[Gac])
				F_ieac = FABints.matrix[Gac][ie][ac];
 
			      value_c += t_jkbe * F_ieac;
			    }

			    /* +t_JkEc * F_IEAB */
			    Ge = Gj ^ Gk ^ Gc;
			    for(e=0; e < avirtpi[Ge]; e++) {
			      E = avir_off[Ge] + e;

			      ec = T2AB.params->colidx[E][C];
			      ie = FAAints.params->rowidx[I][E];

			      t_jkec = F_ieab = 0.0;

			      if(T2AB.params->rowtot[Gjk] && T2AB.params->coltot[Gjk])
				t_jkec = T2AB.matrix[Gjk][jk][ec];

			      if(FAAints.params->rowtot[Gba] && FAAints.params->coltot[Gba])
				F_ieab = FAAints.matrix[Gba][ie][ab];
 
			      value_c += t_jkec * F_ieab;
			    }

			    /* +t_IkAe * F_JeBc */
			    Ge = Gi ^ Gk ^ Ga;
			    for(e=0; e < bvirtpi[Ge]; e++) {
			      E = bvir_off[Ge] + e;

			      ae = T2AB.params->colidx[A][E];
			      je = FABints.params->rowidx[J][E];

			      t_ikae = F_jebc = 0.0;

			      if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
				t_ikae = T2AB.matrix[Gik][ik][ae];

			      if(FABints.params->rowtot[Gbc] && FABints.params->coltot[Gbc])
				F_jebc = FABints.matrix[Gbc][je][bc];
 
			      value_c += t_ikae * F_jebc;
			    }

			    /* -t_IkBe * F_JeAc */
			    Ge = Gi ^ Gk ^ Gb;
			    for(e=0; e < bvirtpi[Ge]; e++) {
			      E = bvir_off[Ge] + e;

			      be = T2AB.params->colidx[B][E];
			      je = FABints.params->rowidx[J][E];

			      t_ikbe = F_jeac = 0.0;

			      if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
				t_ikbe = T2AB.matrix[Gik][ik][be];

			      if(FABints.params->rowtot[Gac] && FABints.params->coltot[Gac])
				F_jeac = FABints.matrix[Gac][je][ac];
 
			      value_c -= t_ikbe * F_jeac;
			    }

			    /* -t_IkEc * F_JEAB */
			    Ge = Gi ^ Gk ^ Gc;
			    for(e=0; e < avirtpi[Ge]; e++) {
			      E = avir_off[Ge] + e;

			      ec = T2AB.params->colidx[E][C];
			      je = FAAints.params->rowidx[J][E];

			      t_ikec = F_jeab = 0.0;

			      if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
				t_ikec = T2AB.matrix[Gik][ik][ec];

			      if(FAAints.params->rowtot[Gba] && FAAints.params->coltot[Gba])
				F_jeab = FAAints.matrix[Gba][je][ab];
 
			      value_c -= t_ikec * F_jeab;
			    }

			    /* +t_IJAE * F_kEcB */
			    Ge = Gi ^ Gj ^ Ga;
			    for(e=0; e < avirtpi[Ge]; e++) {
			      E = avir_off[Ge] + e;

			      ae = T2AA.params->colidx[A][E];
			      ke = FBAints.params->rowidx[K][E];

			      t_ijae = F_kecb = 0.0;

			      if(T2AA.params->rowtot[Gij] && T2AA.params->coltot[Gij])
				t_ijae = T2AA.matrix[Gij][ij][ae];

			      if(FBAints.params->rowtot[Gbc] && FBAints.params->coltot[Gbc])
				F_kecb = FBAints.matrix[Gbc][ke][cb];
 
			      value_c += t_ijae * F_kecb;
			    }

			    /* -t_IJBE * F_kEcA */
			    Ge = Gi ^ Gj ^ Gb;
			    for(e=0; e < avirtpi[Ge]; e++) {
			      E = avir_off[Ge] + e;

			      be = T2AA.params->colidx[B][E];
			      ke = FBAints.params->rowidx[K][E];

			      t_ijbe = F_keca = 0.0;

			      if(T2AA.params->rowtot[Gij] && T2AA.params->coltot[Gij])
				t_ijbe = T2AA.matrix[Gij][ij][be];

			      if(FBAints.params->rowtot[Gac] && FBAints.params->coltot[Gac])
				F_keca = FBAints.matrix[Gac][ke][ca];
 
			      value_c -= t_ijbe * F_keca;
			    }

			    /** <oo||ov> --> connected triples **/

			    /* +t_ImBc * E_kJmA */
			    Gm = Gi ^ Gb ^ Gc;
			    for(m=0; m < boccpi[Gm]; m++) {
			      M = bocc_off[Gm] + m;

			      im = T2AB.params->rowidx[I][M];
			      ma = EBAints.params->colidx[M][A];

			      t_imbc = E_kjma = 0.0;

			      if(T2AB.params->rowtot[Gbc] && T2AB.params->coltot[Gbc])
				t_imbc = T2AB.matrix[Gbc][im][bc];

			      if(EBAints.params->rowtot[Gjk] && EBAints.params->coltot[Gjk])
				E_kjma = EBAints.matrix[Gjk][kj][ma];

			      value_c += t_imbc * E_kjma;
			    }

			    /* -t_ImAc * E_kJmB */
			    Gm = Gi ^ Ga ^ Gc;
			    for(m=0; m < boccpi[Gm]; m++) {
			      M = bocc_off[Gm] + m;

			      im = T2AB.params->rowidx[I][M];
			      mb = EBAints.params->colidx[M][B];

			      t_imac = E_kjmb = 0.0;

			      if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
				t_imac = T2AB.matrix[Gac][im][ac];

			      if(EBAints.params->rowtot[Gjk] && EBAints.params->coltot[Gjk])
				E_kjmb = EBAints.matrix[Gjk][kj][mb];

			      value_c -= t_imac * E_kjmb;
			    }

			    /* +t_IMBA * E_JkMc */
			    Gm = Gi ^ Gb ^ Ga;
			    for(m=0; m < aoccpi[Gm]; m++) {
			      M = aocc_off[Gm] + m;

			      im = T2AA.params->rowidx[I][M];
			      mc = EABints.params->colidx[M][C];

			      t_imba = E_jkmc = 0.0;

			      if(T2AA.params->rowtot[Gba] && T2AA.params->coltot[Gba])
				t_imba = T2AA.matrix[Gba][im][ba];

			      if(EABints.params->rowtot[Gjk] && EABints.params->coltot[Gjk])
				E_jkmc = EABints.matrix[Gjk][jk][mc];

			      value_c += t_imba * E_jkmc;
			    }

			    /* -t_JmBc * E_kImA */
			    Gm = Gj ^ Gb ^ Gc;
			    for(m=0; m < boccpi[Gm]; m++) {
			      M = bocc_off[Gm] + m;

			      jm = T2AB.params->rowidx[J][M];
			      ma = EBAints.params->colidx[M][A];

			      t_jmbc = E_kima = 0.0;

			      if(T2AB.params->rowtot[Gbc] && T2AB.params->coltot[Gbc])
				t_jmbc = T2AB.matrix[Gbc][jm][bc];

			      if(EBAints.params->rowtot[Gik] && EBAints.params->coltot[Gik])
				E_kima = EBAints.matrix[Gik][ki][ma];

			      value_c -= t_jmbc * E_kima;
			    }

			    /* +t_JmAc * E_kImB */
			    Gm = Gj ^ Ga ^ Gc;
			    for(m=0; m < boccpi[Gm]; m++) {
			      M = bocc_off[Gm] + m;

			      jm = T2AB.params->rowidx[J][M];
			      mb = EBAints.params->colidx[M][B];

			      t_jmac = E_kimb = 0.0;

			      if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
				t_jmac = T2AB.matrix[Gac][jm][ac];

			      if(EBAints.params->rowtot[Gik] && EBAints.params->coltot[Gik])
				E_kimb = EBAints.matrix[Gik][ki][mb];

			      value_c += t_jmac * E_kimb;
			    }

			    /* -t_JMBA * E_IkMc */
			    Gm = Gj ^ Gb ^ Ga;
			    for(m=0; m < aoccpi[Gm]; m++) {
			      M = aocc_off[Gm] + m;

			      jm = T2AA.params->rowidx[J][M];
			      mc = EABints.params->colidx[M][C];

			      t_jmba = E_ikmc = 0.0;

			      if(T2AA.params->rowtot[Gba] && T2AA.params->coltot[Gba])
				t_jmba = T2AA.matrix[Gba][jm][ba];

			      if(EABints.params->rowtot[Gik] && EABints.params->coltot[Gik])
				E_ikmc = EABints.matrix[Gik][ik][mc];

			      value_c -= t_jmba * E_ikmc;
			    }

			    /* -t_MkBc * E_JIMA */
			    Gm = Gk ^ Gb ^ Gc;
			    for(m=0; m < aoccpi[Gm]; m++) {
			      M = aocc_off[Gm] + m;

			      mk = T2AB.params->rowidx[M][K];
			      ma = EAAints.params->colidx[M][A];

			      t_mkbc = E_jima = 0.0;

			      if(T2AB.params->rowtot[Gbc] && T2AB.params->coltot[Gbc])
				t_mkbc = T2AB.matrix[Gbc][mk][bc];

			      if(EAAints.params->rowtot[Gji] && EAAints.params->coltot[Gji])
				E_jima = EAAints.matrix[Gji][ji][ma];

			      value_c -= t_mkbc * E_jima;
			    }

			    /* +t_MkAc * E_JIMB */
			    Gm = Gk ^ Ga ^ Gc;
			    for(m=0; m < aoccpi[Gm]; m++) {
			      M = aocc_off[Gm] + m;

			      mk = T2AB.params->rowidx[M][K];
			      mb = EAAints.params->colidx[M][B];

			      t_mkac = E_jimb = 0.0;

			      if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
				t_mkac = T2AB.matrix[Gac][mk][ac];

			      if(EAAints.params->rowtot[Gji] && EAAints.params->coltot[Gji])
				E_jimb = EAAints.matrix[Gji][ji][mb];

			      value_c += t_mkac * E_jimb;
			    }

			    /** disconnected triples **/

			    value_d = 0.0;

			    /* +t_IA * D_JkBc */
			    if(Gi == Ga && Gjk == Gbc) {
			      t_ia = D_jkbc = 0.0;

			      if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
				t_ia = T1A.matrix[Gi][i][a];

			      if(DABints.params->rowtot[Gjk] && DABints.params->coltot[Gjk])
				D_jkbc = DABints.matrix[Gjk][jk][bc];

			      value_d += t_ia * D_jkbc;
			    }

			    /* -t_IB * D_JkAc */
			    if(Gi == Gb && Gjk == Gac) {
			      t_ib = D_jkac = 0.0;

			      if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
				t_ib = T1A.matrix[Gi][i][b];

			      if(DABints.params->rowtot[Gjk] && DABints.params->coltot[Gjk])
				D_jkac = DABints.matrix[Gjk][jk][ac];

			      value_d -= t_ib * D_jkac;
			    }

			    /* -t_JA * D_IkBc */
			    if(Gj == Ga && Gik == Gbc) {
			      t_ja = D_ikbc = 0.0;

			      if(T1A.params->rowtot[Gj] && T1A.params->coltot[Gj])
				t_ja = T1A.matrix[Gj][j][a];

			      if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
				D_ikbc = DABints.matrix[Gik][ik][bc];

			      value_d -= t_ja * D_ikbc;
			    }

			    /* +t_JB * D_IkAc */
			    if(Gj == Gb && Gik == Gac) {
			      t_jb = D_ikac = 0.0;

			      if(T1A.params->rowtot[Gj] && T1A.params->coltot[Gj])
				t_jb = T1A.matrix[Gj][j][b];

			      if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
				D_ikac = DABints.matrix[Gik][ik][ac];

			      value_d += t_jb * D_ikac;
			    }

			    /* +t_kc * D_JIBA */
			    if(Gk == Gc && Gji == Gba) {
			      t_kc = D_jiba = 0.0;

			      if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
				t_kc = T1B.matrix[Gk][k][c];

			      if(DAAints.params->rowtot[Gji] && DAAints.params->coltot[Gji])
				D_jiba = DAAints.matrix[Gji][ji][ba];

			      value_d += t_kc * D_jiba;
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
			    if(fIJ.params->rowtot[Gj])
			      denom += fIJ.matrix[Gj][j][j];
			    if(fij.params->rowtot[Gk])
			      denom += fij.matrix[Gk][k][k];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];
			    if(fab.params->rowtot[Gc])
			      denom -= fab.matrix[Gc][c][c];

			    /*
			    if(fabs(value_c) > 1e-7)
			      fprintf(outfile, "%d %d %d %d %d %d %20.15f\n", I,J,K,A,B,C,value_c);
			    */

			    ET_AAB += (value_d + value_c) * value_c / denom;

			    /* } */ /* A >= B */

			} /* c */
		      } /* b */
		    } /* a */

		  } /* Gb */
		} /* Ga */
	      
	      } /* I >= J */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  ET_AAB /= 2;

  /*  fprintf(outfile, "cnt = %d\n", cnt); */
  /* fprintf(outfile, "ET_AAB = %20.14f\n", ET_AAB); */

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2AA, h);
    dpd_buf4_mat_irrep_close(&T2AB, h);
    dpd_buf4_mat_irrep_close(&FAAints, h);
    dpd_buf4_mat_irrep_close(&FABints, h);
    dpd_buf4_mat_irrep_close(&FBAints, h);
    dpd_buf4_mat_irrep_close(&EAAints, h);
    dpd_buf4_mat_irrep_close(&EABints, h);
    dpd_buf4_mat_irrep_close(&EBAints, h);
    dpd_buf4_mat_irrep_close(&DAAints, h);
    dpd_buf4_mat_irrep_close(&DABints, h);
  }

  dpd_buf4_close(&T2AA);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&FAAints);
  dpd_buf4_close(&FABints);
  dpd_buf4_close(&FBAints);
  dpd_buf4_close(&EAAints);
  dpd_buf4_close(&EABints);
  dpd_buf4_close(&EBAints);
  dpd_buf4_close(&DAAints);
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

  return ET_AAB;
}
