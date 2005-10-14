#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void cc2_sigmaSS(int i, int C_irr);

void cc2_sigma(int i, int C_irr) 
{
  dpdfile2 SIA;
  dpdbuf4 SIjAb;
  dpdfile2 CME;
  dpdbuf4 CMnEf;
  dpdfile2 FAE;
  dpdfile2 FMI;
  dpdfile2 FME;
  dpdbuf4 WMbEj;
  dpdbuf4 WAmEf;
  dpdbuf4 WMnIe;
  dpdbuf4 WAbEi;
  dpdbuf4 WMbIj;
  dpdbuf4 Z;
  dpdbuf4 Z2;
  char lbl[32];
  char CME_lbl[32];
  char CMnEf_lbl[32];
  char SIjAb_lbl[32];

  if (params.eom_ref == 0) { /* RHF */

    cc2_sigmaSS(i,C_irr);

    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    dpd_dot24(&FME,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_file2_close(&FME);
    dpd_file2_close(&SIA);

    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_contract442(&CMnEf, &WAmEf, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&CMnEf);
    dpd_file2_close(&SIA);

    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
    dpd_contract442(&WMnIe, &CMnEf, &SIA, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WMnIe);
    dpd_file2_close(&SIA);

    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmaijDS Z(Ij,Ab)");
    dpd_buf4_init(&WMbIj, CC2_HET1, H_IRR, 10, 0, 10, 0, 0, "CC2 WMbIj");
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_contract244(&CME, &WMbIj, &Z, 0, 0, 1, 1.0, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WMbIj);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WmaijDS Z(jI,bA)");
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb,  -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmaijDS Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb,  -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);

    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(Ij,Ab)");
    dpd_buf4_init(&WAbEi, CC2_HET1, H_IRR, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_contract244(&CME, &WAbEi, &Z, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WabejDS Z(jI,bA)");
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);

    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fbe Z(Ij,Ab)");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "fAB");
    dpd_contract424(&CMnEf, &FAE, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);
    dpd_buf4_close(&CMnEf);

    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "FDD_Fbe Z(jI,bA)");
    dpd_buf4_init(&Z2, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fbe Z(jI,bA)");

    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_axpy(&Z2, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&SIjAb);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fmj Z(Ij,Ab)");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "fIJ");
    dpd_contract244(&FMI, &CMnEf, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&FMI);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "FDD_Fmj Z(jI,bA)"); 
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fmj Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);
  }
  else if (params.eom_ref == 1) { /* ROHF */
    printf("ROHF CC2-LR is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
  else { /* UHF */
    printf("UHF CC2-LR is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
}

void cc2_sigmaSS(int i, int C_irr) 
{
  dpdfile2 SIA;
  dpdfile2 CME;
  dpdfile2 FAE;
  dpdfile2 FMI;
  dpdbuf4 WMbEj;
  char lbl[32];

  if (params.eom_ref == 0) { /* RHF */

    /* sigmaSS */

    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);

    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);

    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_file2_close(&FMI);

    /*
    dpd_buf4_init(&WMbEj, CC2_HET1, H_IRR, 10, 10, 10, 10, 0, "CC2 2 W(jb,ME) + W(Jb,Me)");
    */
    dpd_buf4_init(&WMbEj, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "2 W(jb,ME) + W(Jb,Me)");
    dpd_contract422(&WMbEj, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WMbEj);

    dpd_file2_close(&CME);
    dpd_file2_close(&SIA);
  }
  else if (params.eom_ref == 1) { /* ROHF */
    printf("ROHF CC2-LR is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
  else { /* UHF */
    printf("UHF CC2-LR is not currently implemented\n");
    exit(PSI_RETURN_FAILURE);
  }
}
