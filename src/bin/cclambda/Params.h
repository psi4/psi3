/* Input parameters for cclambda */
struct Params {
  int maxiter;
  double convergence;
  int restart;
  long int memory;
  int cachelev;
  int aobasis;
  int ref;
  int ground;
  double L0; /* 1 for ground states, 0 for excited states */
  double R0; /* only matters for excited states - read from CC_INFO */
  double cceom_energy; /* only matters for excited states - read from CC_INFO */
};

