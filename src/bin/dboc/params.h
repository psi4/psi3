
#ifndef _psi3_bin_dboc_params_h_
#define _psi3_bin_dboc_params_h_

typedef struct {

  enum RefType { rhf=1, rohf=2, uhf=3};

  char *wfn;
  RefType reftype;
  double delta;
  int disp_per_coord;

} Params_t;

#endif
