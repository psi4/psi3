#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void rhf_opdm(void);
void uhf_opdm(void);

void opdm(void)
{
  if(params.ref == 0) return(rhf_opdm());
  else if(params.ref == 2) return(uhf_opdm());
}

void rhf_opdm(void)
{
  int h, i, j, a, b;
  int I, J, A, B;
  double trace=0.0;
  dpdfile2 DIJ, DAB;
  dpdbuf4 T2A, T2B;
  double **OPDM;
  
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
  
  fprintf(outfile, "\n\tTrace of OPDM(2)        = %20.15f\n", fabs(trace));

  OPDM = block_matrix(mo.nmo,mo.nmo);  

  dpd_file2_init(&DIJ, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&DIJ);
  dpd_file2_mat_rd(&DIJ);
 
  for(h=0; h<mo.nirreps; h++) {
    for(i=0; i<mo.occpi[h]; i++) {
      I = mo.qt_occ[mo.occ_off[h] + i];
      for(j=0; j<mo.occpi[h]; j++) {
        J = mo.qt_occ[mo.occ_off[h] + j];
        OPDM[I][J] += DIJ.matrix[h][i][j];
        /* Add OPDM(HF) to OPDM(2) */
	if(I==J) {
	  OPDM[I][J] += 2;
	  trace += OPDM[I][J];
	}
      }
    }
  }
   
  dpd_file2_mat_close(&DIJ);
  dpd_file2_close(&DIJ);

  dpd_file2_init(&DAB, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&DAB);
  dpd_file2_mat_rd(&DAB);
    
  for(h=0; h<mo.nirreps; h++) {
    for(a=0; a<mo.virpi[h]; a++) {
      A = mo.qt_vir[mo.vir_off[h] + a];
      for(b=0; b<mo.virpi[h]; b++) {
        B = mo.qt_vir[mo.vir_off[h] + b];
	OPDM[A][B] += DAB.matrix[h][a][b];
	if(A==B) {
	  trace += OPDM[A][B];
	}
      }
    }
  }

  dpd_file2_mat_close(&DAB);
  fprintf(outfile, "\tTrace of MP2 OPDM       = %20.15f\n", trace);
  dpd_file2_close(&DAB);

  if(params.opdm_print || params.print > 5) {
    fprintf(outfile, "\n\tMP2 OPDM (MO)\n");
    print_mat(OPDM,mo.nmo,mo.nmo,outfile);
  }

  if(params.opdm_write) {
    psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_write_entry(PSIF_MO_OPDM,"MO-basis OPDM",(char*)OPDM[0],sizeof(double)*mo.nmo*mo.nmo);
    psio_close(PSIF_MO_OPDM, 1);
  }
 
  free_block(OPDM);
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
  double **AOPDM, **BOPDM;

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
   
  fprintf(outfile,"\n");
  fprintf(outfile,"\tTrace of Alpha OPDM(2)  = %20.15f\n", fabs(traceA));
  fprintf(outfile,"\tTrace of Beta OPDM(2)   = %20.15f\n", fabs(traceB));
    
  AOPDM = block_matrix(mo.nmo,mo.nmo);  

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
   
  for(h=0; h<mo.nirreps; h++) {
    for(i=0; i<mo.aoccpi[h]; i++) {
      I = mo.qt_aocc[mo.aocc_off[h] + i];
      for(j=0; j<mo.aoccpi[h]; j++) {
	J = mo.qt_aocc[mo.aocc_off[h] + j];
	AOPDM[I][J] += D.matrix[h][i][j];
	/* Add AOPDM(HF) to AOPDM(2) */
	if(I==J) {
	  AOPDM[I][J] += 1;
	  traceA += AOPDM[I][J];
	}
      }
    }
  }
   
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  BOPDM = block_matrix(mo.nmo,mo.nmo);  

  dpd_file2_init(&D, CC_OEI, 0, 2, 2, "Dij");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
   
  for(h=0; h<mo.nirreps; h++) {
    for(i=0; i<mo.boccpi[h]; i++) {
      I = mo.qt_bocc[mo.bocc_off[h] + i];
      for(j=0; j<mo.boccpi[h]; j++) {
	J = mo.qt_bocc[mo.bocc_off[h] + j];
	BOPDM[I][J] += D.matrix[h][i][j];
	/* Add BOPDM(HF) to BOPDM(2) */
	if(I==J) {
	  BOPDM[I][J] += 1;
	  traceB += BOPDM[I][J];
	}
      }
    }
  }
   
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);
    
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
    
  for(h=0; h<mo.nirreps; h++) {
    for(a=0; a<mo.avirpi[h]; a++) {
      A = mo.qt_avir[mo.avir_off[h] + a];
      for(b=0; b<mo.avirpi[h]; b++) {
	B = mo.qt_avir[mo.avir_off[h] + b];
	AOPDM[A][B] += D.matrix[h][a][b];
	if(A==B) {
	  traceA += AOPDM[A][B];
	}
      }
    }
  }

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 3, 3, "Dab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
    
  for(h=0; h<mo.nirreps; h++) {
    for(a=0; a<mo.bvirpi[h]; a++) {
      A = mo.qt_bvir[mo.bvir_off[h] + a];
      for(b=0; b<mo.bvirpi[h]; b++) {
	B = mo.qt_bvir[mo.bvir_off[h] + b];
	BOPDM[A][B] += D.matrix[h][a][b];
	if(A==B) traceB += BOPDM[A][B];
      }
    }
  }
    
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);
    
  fprintf(outfile, "\tTrace of MP2 Alpha OPDM = %20.15f\n", traceA);
  fprintf(outfile, "\tTrace of MP2 Beta  OPDM = %20.15f\n", traceB);
    
  if(params.opdm_print || params.print > 5) {
    fprintf(outfile, "\n\tMP2 Alpha OPDM (MO)\n");
    print_mat(AOPDM,mo.nmo,mo.nmo,outfile);
    fprintf(outfile, "\n\tMP2 Beta OPDM (MO)\n");
    print_mat(BOPDM,mo.nmo,mo.nmo,outfile);
  }

  if(params.opdm_write) {
    psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_write_entry(PSIF_MO_OPDM,"MO-basis Alpha OPDM",(char*)AOPDM[0],sizeof(double)*mo.nmo*mo.nmo);
    psio_write_entry(PSIF_MO_OPDM,"MO-basis Beta OPDM",(char*)BOPDM[0],sizeof(double)*mo.nmo*mo.nmo);
    psio_close(PSIF_MO_OPDM, 1);
  }

  free_block(AOPDM);
  free_block(BOPDM);
}
