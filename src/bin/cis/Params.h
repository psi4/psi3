/*! \file Params.h
    \ingroup (CIS)
    \brief Enter brief description of file here 
*/
/* Input parameters */
struct Params {
  long int memory;
  char *wfn;
  char *diag_method;
  double convergence;
  int maxiter;
  int ref;
  int cis_ref;
  int print;
  int *rpi;
  int local;
};

