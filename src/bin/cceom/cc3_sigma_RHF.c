/*! \file 
    \ingroup (CCEOM)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

/* 
  This function computes contributions to singles and doubles of
  matrix elements of triples:
    SIA   <-- <S|(Dints)           <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt) 
    SIjAb <-- <D|(FME,WAmEf,WMnIe) <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt) 
  Irrep variables are:
    Sirr <--       Wirr                  WX3irr  ^ Cirr
                                  (           X3irr            )
  These are used to make X3 quantity in T3_RHF:
    CIjAb, WAbEi, WMbIj, fIJ2, fAB2, omega
*/

void cc3_sigma_RHF_obsolete(dpdbuf4 *CIjAb, dpdbuf4 *WAbEi, dpdbuf4 *WMbIj,
    int do_singles, dpdbuf4 *Dints, dpdfile2 *SIA,
    int do_doubles, dpdfile2 *FME, dpdbuf4 *WAmEf, dpdbuf4 *WMnIe, 
    dpdbuf4 *SIjAb, double energy)
{
  int h, nirreps;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int Gi, Gj, Gk, Gl, Ga, Gb, Gc, Gd;
  int i, j, k, l, a, b, c, d;
  int I, J, K, L, A, B, C, D;
  int kj, jk, ji, ij, ik, ki;
  int Gkj, Gjk, Gji, Gij, Gik, Gki;
  int Gijk, Sirr, Cirr, WX3irr, Wirr, X3irr;
  int ab, ba, ac, ca, bc, cb;
  int Gab, Gba, Gac, Gca, Gbc, Gcb;
  int id, jd, kd, ad, bd, cd;
  int il, jl, kl, la, lb, lc, li, lk;
  int da, di, dj, dk;
  int Gad, Gdi, Gdj, Gdk, Glc, Gli, Glk;
  double value, F_val, t_val, E_val;
  double dijk, denom;
  double value_ia, value_ka, denom_ia, denom_ka;
  dpdfile2 fIJ, fIJ2, fAB, fAB2, SIA_inc;
  dpdbuf4 SIjAb_inc;
  double ***T3, ***W3;
  int nv;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  /* these are sent to T3 function */
  dpd_file2_init(&fIJ2, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB2, CC_OEI, 0, 1, 1, "fAB");

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_init(FME);
  dpd_file2_mat_rd(FME);

  Cirr = CIjAb->file.my_irrep;
  WX3irr = WAbEi->file.my_irrep;
  X3irr = WX3irr^Cirr;
  Wirr = WAmEf->file.my_irrep;
  Sirr = SIjAb->file.my_irrep;
  if (Sirr != (X3irr^Wirr)) {
    fprintf(outfile,"problem with irreps in cc3_sigma_RHF()\n"); 
    exit(1);
  }
  dpd_file2_init(&SIA_inc, CC_TMP0, Sirr, 0, 1, "CC3 SIA");  /* T3->S1 increment */
  dpd_file2_mat_init(&SIA_inc);
  dpd_buf4_init(&SIjAb_inc, CC_TMP0, Sirr, 0, 5, 0, 5, 0, "CC3 SIjAb");
  dpd_buf4_scm(&SIjAb_inc, 0);

  for(h=0; h < nirreps; h++) {
    if (do_singles) {
      dpd_buf4_mat_irrep_init(Dints, h);
      dpd_buf4_mat_irrep_rd(Dints, h);
    }
    dpd_buf4_mat_irrep_init(WAmEf, h);
    dpd_buf4_mat_irrep_rd(WAmEf, h);
    dpd_buf4_mat_irrep_init(WMnIe, h);
    dpd_buf4_mat_irrep_rd(WMnIe, h);
    dpd_buf4_mat_irrep_init(&SIjAb_inc, h);
  }

  for(h=0,nv=0; h < nirreps; h++) nv += virtpi[h];

  W3 = (double ***) malloc(nirreps * sizeof(double **));
  T3 = init_3d_array(nv, nv, nv);

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gji = Gi ^ Gj; 
      for(Gk=0; Gk < nirreps; Gk++) {
        Gkj = Gjk = Gk ^ Gj;
        Gik = Gki = Gi ^ Gk;
        Gijk = Gi ^ Gj ^ Gk;
        
        /* allocate memory for all irrep blocks of (ab,c) */
        for(Gab=0; Gab < nirreps; Gab++) {
          Gc = Gab ^ Gijk ^ X3irr;
          W3[Gab] = dpd_block_matrix(WAbEi->params->coltot[Gab], virtpi[Gc]);
        }
        
        for(i=0; i < occpi[Gi]; i++) {
          I = occ_off[Gi] + i;
          for(j=0; j < occpi[Gj]; j++) {
            J = occ_off[Gj] + j;
            for(k=0; k < occpi[Gk]; k++) {
              K = occ_off[Gk] + k;
              
              ij = CIjAb->params->rowidx[I][J];
              ji = CIjAb->params->rowidx[J][I];
              ik = CIjAb->params->rowidx[I][K];
              ki = CIjAb->params->rowidx[K][I];
              jk = CIjAb->params->rowidx[J][K];
              kj = CIjAb->params->rowidx[K][J];
              
              T3_RHF(W3, nirreps, I, Gi, J, Gj, K, Gk, CIjAb, WAbEi, WMbIj, 
                     &fIJ2, &fAB2, occpi, occ_off, virtpi, vir_off, energy);

              /* sort (ab,c) into T3[A][B][C] */
              for(Gab=0; Gab < nirreps; Gab++) {
                Gc = Gab ^ Gijk ^ X3irr;
                for (ab=0; ab < WAbEi->params->coltot[Gab]; ++ab) {
                  A = WAbEi->params->colorb[Gab][ab][0];
                  B = WAbEi->params->colorb[Gab][ab][1];
                  Ga = WAbEi->params->rsym[A];
                  Gb = Ga ^ Gab;
                  /* a = A - vir_off[Ga];
                  b = B - vir_off[Gb]; */ 
                  for (c=0; c<virtpi[Gc]; ++c) {
                    C = c + vir_off[Gc];
                    T3[A][B][C] = W3[Gab][ab][c];
                  }
                }
              }
              /*** X3 --> SIA Contributions ***/
              if (do_singles) {
                for(Ga=0; Ga < nirreps; Ga++) {
                  for(a=0; a < virtpi[Ga]; a++) {
                    A = vir_off[Ga] + a;
  
                    value_ia = 0.0;
                    value_ka = 0.0;
                    for(Gb=0; Gb < nirreps; Gb++) {
                      for(b=0; b < virtpi[Gb]; b++) {
                        B = vir_off[Gb] + b;
  
                        Gc = Gijk ^ Ga ^ Gb ^ X3irr;
                        Gbc = Gb ^ Gc;
  
                        for(c=0; c < virtpi[Gc]; c++) {
                          C = vir_off[Gc] + c;
  
                          bc = Dints->params->colidx[B][C];
  
                          if(Gi^Ga==Sirr && Gjk == Gbc) {
  
                            if(Dints->params->rowtot[Gjk] && Dints->params->coltot[Gjk])
                              value_ia += T3[A][B][C] * Dints->matrix[Gjk][jk][bc];
                          }
  
                          if(Gk^Ga==Sirr && Gji == Gbc) {
  
                            if(Dints->params->rowtot[Gji] && Dints->params->coltot[Gji])
                              value_ka -= T3[A][B][C] * Dints->matrix[Gji][ji][bc];
                          }
  
                        } /* c */
                      } /* b */
                    } /* Gb */
  
                    if(SIA_inc.params->rowtot[Gi] && SIA_inc.params->coltot[Gi^Sirr])
                      SIA_inc.matrix[Gi][i][a] += value_ia;
                    if(SIA_inc.params->rowtot[Gk] && SIA_inc.params->coltot[Gk^Sirr])
                      SIA_inc.matrix[Gk][k][a] += value_ka;
  
                  } /* a */
                } /* Ga */
              } /*** end X3 --> SIA ***/

              /*** X3 --> SIjAb Contributions ***/
              if (do_doubles) {
                for(Ga=0; Ga < nirreps; Ga++) {
                  for(Gb=0; Gb < nirreps; Gb++) {
                    Gab = Ga ^ Gb;
                    Gc = Gijk ^ Gab ^ X3irr;
  
                    for(a=0; a < virtpi[Ga]; a++) {
                      A = vir_off[Ga] + a;
                      for(b=0; b < virtpi[Gb]; b++) {
                        B = vir_off[Gb] + b;
  
                        ab = SIjAb_inc.params->colidx[A][B];
                        ba = SIjAb_inc.params->colidx[B][A];
  
                        if( ((Gij^Gab)==Sirr) && ((Gk^Gc)==Wirr) ) {
  
                          value = 0.0;
                          for(c=0; c < virtpi[Gc]; c++) {
                            C = vir_off[Gc] + c;
                            if(FME->params->rowtot[Gk] && FME->params->coltot[Gc^Wirr])
                              value += T3[A][B][C] * FME->matrix[Gk][k][c];
                          }
  
                          if(SIjAb_inc.params->rowtot[Gij] && SIjAb_inc.params->coltot[Gab^Sirr]) {
                            SIjAb_inc.matrix[Gij][ij][ab] += value;
                            SIjAb_inc.matrix[Gij][ji][ba] += value;
                          }
                        }
                        if( ((Gjk^Gab)==Sirr) && ((Gi^Gc)==Wirr) ) {
  
                          value = 0.0;
                          for(c=0; c < virtpi[Gc]; c++) {
                            C = vir_off[Gc] + c;
                            if(FME->params->rowtot[Gi] && FME->params->coltot[Gc^Wirr])
                              value -= T3[A][B][C] * FME->matrix[Gi][i][c];
                          }
  
                          if(SIjAb_inc.params->rowtot[Gjk] && SIjAb_inc.params->coltot[Gab^Sirr]) {
                            SIjAb_inc.matrix[Gjk][kj][ab] += value;
                            SIjAb_inc.matrix[Gjk][jk][ba] += value;
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
  
                      di = WAmEf->params->rowidx[D][I];
                      dj = WAmEf->params->rowidx[D][J];
                      dk = WAmEf->params->rowidx[D][K];
  
                      for(a=0; a < virtpi[Ga]; a++) {
                        A = vir_off[Ga] + a;
  
                        ad = SIjAb_inc.params->colidx[A][D];
                        da = SIjAb_inc.params->colidx[D][A];
  
                        if( (Gij^Gad)==Sirr) {
                          value = 0.0;
                          for(Gb=0; Gb < nirreps; Gb++) {
                            Gc = Gijk ^ Ga ^ Gb ^ X3irr;
                            Gbc = Gb ^ Gc;
                            if( (Gdk^Gbc)==Wirr) {
                              for(b=0; b < virtpi[Gb]; b++) {
                                B = vir_off[Gb] + b;
                                for(c=0; c < virtpi[Gc]; c++) {
                                  C = vir_off[Gc] + c;
  
                                  bc = WAmEf->params->colidx[B][C];
  
                                  if(WAmEf->params->rowtot[Gdk] && WAmEf->params->coltot[Gdk^Wirr])
                                    value += 2.0 * T3[A][B][C] * WAmEf->matrix[Gdk][dk][bc];
                                }
                              }
                            }
                          }
  
                          if(SIjAb_inc.params->rowtot[Gij] && SIjAb_inc.params->coltot[Gij^Sirr]) {
                            SIjAb_inc.matrix[Gij][ij][ad] += value;
                            SIjAb_inc.matrix[Gij][ji][da] += value;
                          }
                        }
                        if((Gjk^Gad)==Sirr) {
                          value = 0.0;
                          for(Gb=0; Gb < nirreps; Gb++) {
                            Gc = Gijk ^ Ga ^ Gb ^ X3irr;
                            Gbc = Gb ^ Gc;
                            if((Gdi^Gbc)==Wirr) {
                              for(b=0; b < virtpi[Gb]; b++) {
                                B = vir_off[Gb] + b;
                                for(c=0; c < virtpi[Gc]; c++) {
                                  C = vir_off[Gc] + c;
  
                                  bc = WAmEf->params->colidx[B][C];
  
                                  if(WAmEf->params->rowtot[Gdi] && WAmEf->params->coltot[Gdi^Wirr])
                                    value -= T3[A][B][C] * WAmEf->matrix[Gdi][di][bc];
                                }
                              }
                            }
                          }
  
                          if(SIjAb_inc.params->rowtot[Gjk] && SIjAb_inc.params->coltot[Gjk^Sirr]) {
                            SIjAb_inc.matrix[Gjk][kj][ad] += value;
                            SIjAb_inc.matrix[Gjk][jk][da] += value;
                          }
  
                        }
                        if((Gik^Gad)==Sirr) {
                          value = 0.0;
                          for(Gb=0; Gb < nirreps; Gb++) {
                            Gc = Gijk ^ Ga ^ Gb ^ X3irr;
                            Gbc = Gb ^ Gc;
                            if((Gdj^Gbc)==Wirr) {
                              for(b=0; b < virtpi[Gb]; b++) {
                                B = vir_off[Gb] + b;
                                for(c=0; c < virtpi[Gc]; c++) {
                                  C = vir_off[Gc] + c;
  
                                  bc = WAmEf->params->colidx[B][C];
  
                                  if(WAmEf->params->rowtot[Gdj] && WAmEf->params->coltot[Gdj^Wirr])
                                    value -= T3[A][B][C] * WAmEf->matrix[Gdj][dj][bc];
                                }
                              }
                            }
                          }
  
                          if(SIjAb_inc.params->rowtot[Gik] && SIjAb_inc.params->coltot[Gik^Sirr]) {
                            SIjAb_inc.matrix[Gki][ki][da] += value;
                            SIjAb_inc.matrix[Gki][ik][ad] += value;
                          }
                        }
                      } 
                    } 
                  } 
                } /* end Wamef terms */
                for(Gl=0; Gl < nirreps; Gl++) {
                  Gli = Gl ^ Gi;
                  Glk = Gl ^ Gk;
  
                  for(Ga=0; Ga < nirreps; Ga++) {
                    for(Gb=0; Gb < nirreps; Gb++) {
                      Gab = Ga ^ Gb;
  
                      if((Gli^Gab)==Sirr) {
  
                        for(l=0; l < occpi[Gl]; l++) {
                          L = occ_off[Gl] + l;
  
                          li = SIjAb_inc.params->rowidx[L][I];
                          il = SIjAb_inc.params->rowidx[I][L];
  
                          for(a=0; a < virtpi[Ga]; a++) {
                            A = vir_off[Ga] + a;
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
  
                              ab = SIjAb_inc.params->colidx[A][B];
                              ba = SIjAb_inc.params->colidx[B][A];
  
                              value = 0.0;
                              for(Gc=0; Gc < nirreps; Gc++) {
                                Glc = Gl ^ Gc;
  
                                if((Gjk^Glc)==Wirr) {
  
                                  for(c=0; c < virtpi[Gc]; c++) {
                                    C = vir_off[Gc] + c;
                                    lc = WMnIe->params->colidx[L][C];
  
                                    if(WMnIe->params->rowtot[Gjk] && WMnIe->params->coltot[Gjk^Wirr])
                                      value -= T3[A][B][C] * WMnIe->matrix[Gjk][jk][lc];
                                  }
                                }
  
                              }
  
                              value *= 2.0;
  
                              if(SIjAb_inc.params->rowtot[Gli] && SIjAb_inc.params->coltot[Gli^Sirr]) {
                                SIjAb_inc.matrix[Gli][li][ba] += value;
                                SIjAb_inc.matrix[Gli][il][ab] += value;
                              }
                            } 
                          } 
                        } 
                      } 
                      if((Glk^Gab)==Sirr) {
  
                        for(l=0; l < occpi[Gl]; l++) {
                          L = occ_off[Gl] + l;
  
                          lk = SIjAb_inc.params->rowidx[L][K];
                          kl = SIjAb_inc.params->rowidx[K][L];
  
                          for(a=0; a < virtpi[Ga]; a++) {
                            A = vir_off[Ga] + a;
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
  
                              ab = SIjAb_inc.params->colidx[A][B];
                              ba = SIjAb_inc.params->colidx[B][A];
  
                              value = 0.0;
                              for(Gc=0; Gc < nirreps; Gc++) {
                                Glc = Gl ^ Gc;
  
                                if((Gji^Glc)==Wirr) {
  
                                  for(c=0; c < virtpi[Gc]; c++) {
                                    C = vir_off[Gc] + c;
                                    lc = WMnIe->params->colidx[L][C];
  
                                    if(WMnIe->params->rowtot[Gji] && WMnIe->params->coltot[Gji^Wirr])
                                      value += T3[A][B][C] * WMnIe->matrix[Gji][ji][lc];
                                  }
                                }
  
                              } 
  
                              if(SIjAb_inc.params->rowtot[Glk] && SIjAb_inc.params->coltot[Glk^Sirr]) {
                                SIjAb_inc.matrix[Glk][lk][ba] += value;
                                SIjAb_inc.matrix[Glk][kl][ab] += value;
                              }
                            } 
                          } 
                        } 
                      } 

                      if((Gli^Gab)==Sirr) {
  
                        for(l=0; l < occpi[Gl]; l++) {
                          L = occ_off[Gl] + l;
  
                          li = SIjAb_inc.params->rowidx[L][I];
                          il = SIjAb_inc.params->rowidx[I][L];
  
                          for(a=0; a < virtpi[Ga]; a++) {
                            A = vir_off[Ga] + a;
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
  
                              ab = SIjAb_inc.params->colidx[A][B];
                              ba = SIjAb_inc.params->colidx[B][A];
  
                              value = 0.0;
                              for(Gc=0; Gc < nirreps; Gc++) {
                                Glc = Gl ^ Gc;
  
                                if((Gjk^Glc)==Wirr) {
  
                                  for(c=0; c < virtpi[Gc]; c++) {
                                    C = vir_off[Gc] + c;
                                    lc = WMnIe->params->colidx[L][C];
  
                                    if(WMnIe->params->rowtot[Gjk] && WMnIe->params->coltot[Gjk^Wirr])
                                      value += T3[A][B][C] * WMnIe->matrix[Gjk][kj][lc];
                                  }
                                }
  
                              }
  
                              if(SIjAb_inc.params->rowtot[Gli] && SIjAb_inc.params->coltot[Gli^Sirr]) {
                                SIjAb_inc.matrix[Gli][li][ba] += value;
                                SIjAb_inc.matrix[Gli][il][ab] += value;
                              }
                            } 
                          } 
                        } 
                      } 
  
                    }
                  }
                } /* end Wmnie terms */
              } /*** end do_doubles */
            } /* k */
          } /* j */
         } /* i */

        for(Gab=0; Gab < nirreps; Gab++) {
          /* This will need to change for non-totally-symmetric cases */
          Gc = Gab ^ Gijk;
          dpd_free_block(W3[Gab], WAbEi->params->coltot[Gab^WX3irr], virtpi[Gc]);
        }
      } /* Gk */
    } /* Gj */
  } /* Gi */
  free_3d_array(T3, nv, nv);

  /* close up files and update sigma vectors */
  for(h=0; h < nirreps; h++) {
    if (do_singles) dpd_buf4_mat_irrep_close(Dints, h);
    dpd_buf4_mat_irrep_close(WAmEf, h);
    dpd_buf4_mat_irrep_close(WMnIe, h);
    dpd_buf4_mat_irrep_wrt(&SIjAb_inc, h);
    dpd_buf4_mat_irrep_close(&SIjAb_inc, h);
  }
  dpd_file2_mat_wrt(&SIA_inc);
  dpd_file2_mat_close(&SIA_inc);
  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fIJ2);
  dpd_file2_close(&fAB2);

  /*
  fprintf(outfile, "triples S1norm = %20.10f\n", sqrt(dpd_file2_dot_self(&SIA_inc)));
  fprintf(outfile, "triples S2norm = %20.10f\n", sqrt(dpd_buf4_dot_self(&SIjAb_inc)));
  */
  if (do_singles) dpd_file2_axpy(&SIA_inc, SIA, 1, 0);
  dpd_file2_close(&SIA_inc);
  dpd_buf4_axpy(&SIjAb_inc, SIjAb, 1);
  dpd_buf4_close(&SIjAb_inc);
}
