struct Params {
  int print;             /* Output level control */
  long int memory;       /* Memory available (in bytes) */
  int cachelev;          /* cacheing level for libdpd */
  int ref;               /* reference determinant (0=RHF, 1=ROHF, 2=UHF) */
  double omega;          /* energy of applied field (a.u) for dynamic polarizabilities */
  int maxiter;           /* maximum number of iterations allowed to converge perturbed amp eqns. */
  double convergence;    /* convergence criterion for perturbed wfns */
  int diis;              /* boolean for using DIIS extrapolation */
  char *prop;            /* user-selected property */
  int local;             /* boolean for simluation of local correlation */
};
