#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void rhf_opdm(void);
void uhf_opdm(void);

void opdm(void)
{
  if(params.ref == 0) rhf_opdm();
  else if(params.ref == 2) uhf_opdm();
}

void rhf_opdm(void)
{
  int h, i, j, a, b;
  int I, J, A, B;
  double trace=0.0;
  dpdfile2 DIJ, DAB;
  dpdbuf4 T2A, T2B;
  
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_scmcopy(&T2A, CC_TAMPS, "2 tIjAb - tIjBa", 2);
  dpd_buf4_sort_axpy(&T2A, CC_TAMPS, pqsr, 0, 5, "2 tIjAb - tIjBa", -1);
  dpd_buf4_close(&T2A);
  
  dpd_file2_init(&DIJ, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_contract442(&T2A, &T2B, &DIJ, 0, 0, -2, 0);
  dpd_buf4_close(&T2B);
  dpd_buf4_close(&T2A);
  trace += dpd_file2_trace(&DIJ);
  dpd_file2_close(&DIJ);
 
  dpd_file2_init(&DAB, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_contract442(&T2A, &T2B, &DAB, 3, 3, 2, 0);
  dpd_buf4_close(&T2B);
  dpd_buf4_close(&T2A);
  trace += dpd_file2_trace(&DAB);
  dpd_file2_close(&DAB);
  
  /*fprintf(outfile, "\n\tTrace of OPDM(2)        = %20.15f\n", fabs(trace));*/

  return;
}

void uhf_opdm(void)
{
  int h, i, j, a, b;
  int I, J, A, B;
  double traceA=0.0;
  double traceB=0.0;
  dpdfile2 D;
  dpdfile2 T1A, T1B;
  dpdbuf4 T2A, T2B;

  if(params.semicanonical) {
    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1A, &T1B, &D, 0, 0, -1, 0);
    dpd_file2_close(&T1A);
    dpd_file2_close(&T1B);
    
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_contract442(&T2A, &T2B, &D, 0, 0, -1, 1);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  else {	
    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB"); 
    dpd_contract442(&T2A, &T2B, &D, 0, 0, -1, 0);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_contract442(&T2A, &T2B, &D, 0, 0, -1, 1);
  traceA += dpd_file2_trace(&D);
  dpd_buf4_close(&T2A);
  dpd_buf4_close(&T2B);
  dpd_file2_close(&D);
    
  if(params.semicanonical) {
    dpd_file2_init(&D, CC_OEI, 0, 2, 2, "Dij");
    dpd_file2_init(&T1A, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&T1A, &T1B, &D, 0, 0, -1, 0);
    dpd_file2_close(&T1A);
    dpd_file2_close(&T1B);
    
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_contract442(&T2A, &T2B, &D, 0, 0, -1, 1);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  else { 
    dpd_file2_init(&D, CC_OEI, 0, 2, 2, "Dij");
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab"); 
    dpd_contract442(&T2A, &T2B, &D, 0, 0, -1, 0);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_contract442(&T2A, &T2B, &D, 1, 1, -1, 1);
  traceB += dpd_file2_trace(&D);
  dpd_buf4_close(&T2A);
  dpd_buf4_close(&T2B);
  dpd_file2_close(&D);
  
  if(params.semicanonical) {
    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1A, &T1B, &D, 1, 1, 1, 0);
    dpd_file2_close(&T1A);
    dpd_file2_close(&T1B);
    
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_contract442(&T2A, &T2B, &D, 3, 3, 1, 1);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  else {
    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB"); 
    dpd_contract442(&T2A, &T2B, &D, 3, 3, 1, 0);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_contract442(&T2A, &T2B, &D, 2, 2, 1, 1);
  traceA += dpd_file2_trace(&D);
  dpd_buf4_close(&T2A);
  dpd_buf4_close(&T2B);
  dpd_file2_close(&D);
    
  if(params.semicanonical) {
    dpd_file2_init(&D, CC_OEI, 0, 3, 3, "Dab");
    dpd_file2_init(&T1A, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&T1A, &T1B, &D, 1, 1, 1, 0);
    dpd_file2_close(&T1A);
    dpd_file2_close(&T1B);
    
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_contract442(&T2A, &T2B, &D, 3, 3, 1, 1);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  else {
    dpd_file2_init(&D, CC_OEI, 0, 3, 3, "Dab");
    dpd_buf4_init(&T2A, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_buf4_init(&T2B, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab"); 
    dpd_contract442(&T2A, &T2B, &D, 3, 3, 1, 0);
    dpd_buf4_close(&T2A);
    dpd_buf4_close(&T2B);
  }
  
  dpd_buf4_init(&T2A, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2B, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_contract442(&T2A, &T2B, &D, 3, 3, 1, 1);
  traceB += dpd_file2_trace(&D);
  dpd_buf4_close(&T2A);
  dpd_buf4_close(&T2B);
  dpd_file2_close(&D);
   
  /*
  fprintf(outfile,"\n");
  fprintf(outfile,"\tTrace of Alpha OPDM(2)  = %20.15f\n", fabs(traceA));
  fprintf(outfile,"\tTrace of Beta OPDM(2)   = %20.15f\n", fabs(traceB));
  */
    
  return;
}
