/*! \file energy.c
    \ingroup (MP2)
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double rhf_energy(void);
double uhf_energy(void);

double energy(void)
{
  if(params.ref == 0) return(rhf_energy());
  else if(params.ref == 2) return(uhf_energy());
}

double rhf_energy(void) 
{
  double E = 0;
  dpdbuf4 tIjAb;
  dpdbuf4 D;
  dpdbuf4 S;
  double s_energy, t_energy, scs_energy;

  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  E = dpd_buf4_dot(&D, &tIjAb);

  if (params.scs == 1) {
    dpd_buf4_init(&S, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    s_energy = dpd_buf4_dot(&S, &tIjAb);
    dpd_buf4_close(&S);
    t_energy = (E - s_energy);
    fprintf(outfile,"\n\tSinglet correlation energy  = %20.15f\n",s_energy);
    fprintf(outfile,"\tTriplet correlation energy  = %20.15f\n",t_energy);
    
    s_energy = params.scs_scale_s * s_energy;
    t_energy = params.scs_scale_t * t_energy;
    scs_energy = s_energy + t_energy;

    fprintf(outfile,"\n\tScale_S correlation energy  = %20.15f\n",s_energy);
    fprintf(outfile,"\tScale_T correlation energy  = %20.15f\n",t_energy);
    fprintf(outfile,"\tSCS-MP2 total energy        = %20.15f\n",mo.Escf + 
      scs_energy);
  }

  dpd_buf4_close(&tIjAb);
  dpd_buf4_close(&D);


  return(E);
}

double uhf_energy(void)
{
  double E1A = 0; double E1B = 0;	
  double E2AA = 0; double E2BB = 0; double E2AB = 0;
  dpdfile2 T1, F;
  dpdbuf4 tIJAB, tijab, tIjAb,  D;
  
  if(params.semicanonical) {
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    E1A = dpd_file2_dot(&F, &T1);
    dpd_file2_close(&F);
    dpd_file2_close(&T1);

    dpd_file2_init(&F, CC_OEI, 0, 2, 3, "fia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    E1B = dpd_file2_dot(&F, &T1);
    dpd_file2_close(&F);
    dpd_file2_close(&T1);
  }
  
  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = dpd_buf4_dot(&D, &tIJAB);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tIJAB);

  dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = dpd_buf4_dot(&D, &tijab);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tijab);
  
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = dpd_buf4_dot(&D, &tIjAb);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tIjAb);
  
  if(params.semicanonical)
    return(E1A+E1B+E2AA+E2BB+E2AB);
  else
    return(E2AA+E2BB+E2AB);
}

