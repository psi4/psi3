/* Input parameters for cclambda */
struct Params {
  int maxiter;
  double convergence;
  int restart;
  long int memory;
  int cachelev;
  int aobasis;
  int ref;
  int ground; /* boolean, =1 implies only ground state calculation */
  double L0; /* 1 for ground states, 0 for excited states */
  double **R0; /* read from CC_INFO if necessary */
  double **cceom_energy; /* read from CC_INFO if necessary */
  int *states_per_irrep; /* determined from input file */
  int *Ls_per_irrep;
  int local;  /* boolean for using simulated local-CC framework */
};

