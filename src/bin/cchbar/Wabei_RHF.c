#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/** Wabei intermediates are stored here as (ei,ab) **/

void Wabei_RHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;
  double value;

  if(params.print & 2) {
    fprintf(outfile, "\tF -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  dpd_buf4_copy(&F, CC_HBAR, "WEiAb");
  dpd_buf4_close(&F);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** - F_ME t_Mi^Ab **/
  if(params.print & 2) {
    fprintf(outfile, "\tFME*T2 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** <Ab|Ef> t_i^f **/
  if(params.print & 2) {
    fprintf(outfile, "\tB*T1 -> Wabei...");
    fflush(outfile);
  }
  /* NB: this term is sorted into WEiAb after the F*T1*T1 contributions below */
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)"); 
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);
  dpd_buf4_close(&Z);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Prepare intermediates for second Wabef contribution to Wabei **/

  /** t_i^f <Am|Ef> --> Z(Am,Ei) **/
  if(params.print & 2) {
    fprintf(outfile, "\tF*T1*T1 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(Am,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ai|bc>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** - t_m^b Z(mA,Ei) --> Z1(Ab,Ei) **/
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(Am,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  /* NB: high memory req. in this contraction b/c of Z transpose */
  dpd_contract424(&Z, &T1, &W, 1, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** <Mb|Ef> t_i^f --> Z(Mb,Ei) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** - T_M^A Z(Mb,Ei) --> Z(Ab,Ei) **/
  /* NB: This term is sorted into WEiAb at the very last term of the Wabei intermediate */
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z1, 0, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z1);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Final term of Wabef contribution to Wabei **/

  /** t_i^f <Mn|Ef> --> Z(Mn,Ei) **/
  if(params.print & 2) {
    fprintf(outfile, "\tD*T1*T2 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(Mn,Ei)");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&D);
  /** Z(Mn,Ei) Tau(Mn,Ab) --> W(Ei,Ab) **/
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  dpd_contract444(&Z, &T, &W, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&T);
  dpd_buf4_close(&Z);
  if(params.print & 2) fprintf(outfile, "done.\n");


  /** <Mn|Ei> Tau(Mn,Ab) --> W(Ei,Ab) **/
  if(params.print & 2) {
    fprintf(outfile, "\tE*T2 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  dpd_contract444(&E, &T, &W, 1, 1, 1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&T);
  dpd_buf4_close(&E);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** -<Mb|Ef> t_Mi^Af - <MA||EF> t_iM^bF + <mA|fE> t_im^bf --> W(Ei,Ab) **/
  if(params.print & 2) {
    fprintf(outfile, "\tF*T2 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_buf4_init(&Z, CC_TMP1, 0, 5, 10, 5, 10, 0, "Z(Eb,iA)");
  dpd_contract444(&F, &T2, &Z, 1, 1, -1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  value = dpd_buf4_dot_self(&Z);
  fprintf(outfile, "\nZ(Eb,iA) * Z(Eb,iA) = %20.12f\n", value);
  fflush(outfile);
  dpd_buf4_sort_axpy(&Z, CC_HBAR, prsq, 11, 5, "WEiAb", 1);  /***** BOTTLENECK *****/
  dpd_buf4_close(&Z);

  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  value = dpd_buf4_dot_self(&W);
  fprintf(outfile, "\nWEiAb * WEiAb = %20.12f\n", value);
  fflush(outfile);
  dpd_buf4_close(&W);

  Wabei_RHF_FT2_a();

  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_init(&Z, CC_TMP2, 0, 5, 10, 5, 10, 0, "Z(EA,ib)");
  dpd_contract444(&F, &T2, &Z, 1, 1, -1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_sort(&F, CC_TMP2, prsq, 10, 5, "F <ia|bc> (ib,ca)");
  dpd_buf4_close(&F);
  dpd_buf4_init(&F, CC_TMP2, 0, 10, 5, 10, 5, 0, "F <ia|bc> (ib,ca)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_buf4_init(&Z, CC_TMP2, 0, 5, 10, 5, 10, 0, "Z(EA,ib)");
  dpd_contract444(&F, &T2, &Z, 1, 1, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_sort_axpy(&Z, CC_HBAR, prqs, 11, 5, "WEiAb", 1);  /***** BOTTLENECK *****/
  dpd_buf4_close(&Z);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Final terms of Wabei **/

  /** t_in^bf  <Mn|Ef> + t_iN^bF <MN||EF> --> Z1_MEib **/
  if(params.print & 2) {
    fprintf(outfile, "\tT1*(C+D*T2) -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,ib)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_contract444(&D, &T2, &Z, 0, 0, -1, 1);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);
  dpd_buf4_sort(&Z, CC_TMP0, psqr, 10, 11, "Z(Mb,Ei)");
  dpd_buf4_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mE,iA)");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 10, "Z(mA,iE)");
  dpd_buf4_close(&Z);

  /** - t_M^A ( <Mb|Ei> + Z(Mb,Ei) ) --> Z1(Ab,Ei) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  dpd_buf4_axpy(&D, &Z, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z1, 0, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  /* NB: This sort_axpy implicitly includes terms from B*T1 and F*T1*T1 from above */
  dpd_buf4_sort_axpy(&Z1, CC_HBAR, rspq, 11, 5, "WEiAb", 1);   /******** BOTTLENECK!!! *********/
  dpd_buf4_close(&Z1);

  /** t_m^b ( - <mA|iE> + Z(mA,iE) ) --> Z2(Ab,Ei) **/
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mA,iE)");
  dpd_buf4_axpy(&C, &Z, -1.0);
  dpd_buf4_close(&C);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&Z2, CC_TMP2, 0, 5, 10, 5, 10, 0, "Z2(bA,iE)");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, CC_HBAR, srqp, 11, 5, "WEiAb", 1);  /******** BOTTLENECK!!! **********/
  dpd_buf4_close(&Z2);
  if(params.print & 2) fprintf(outfile, "done.\n");

}
