/* Input parameters */
struct Params {
  int maxiter;
  double convergence;
  int restart;
  long int memory;
  char *aobasis;
  int cachelev;
  int cachetype;
  int ref;
  int diis;
  char *wfn;
  int print;
  int local;
  int num_amps;
  int print_mp2_amps;
  int brueckner;
  double bconv;
  int analyze;
  int print_pair_energies;
  int spinadapt_energies;
  int semicanonical;
  int local_mos;
  int dertype;
  int t2_coupled;
};
