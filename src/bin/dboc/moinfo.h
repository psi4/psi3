
#ifndef _psi3_bin_dboc_moinfo_h_
#define _psi3_bin_dboc_moinfo_h_

typedef struct {
  int num_so;
  int num_mo;
  int ndocc;
  int nsocc;
  int nalpha;
  int nbeta;
  int nfzc;
  int nfzv;
  int nact;
  short int* QTS_to_pitzer;
} MOInfo_t;

#endif

