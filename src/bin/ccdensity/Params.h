/* Input parameters for cclambda */
struct Params {
  double tolerance;
  long int memory;
  int cachelev;
  int aobasis;
  int ref;
  int relax_opdm;
  int ground;
  int user_transition; /* was L specified on command-line? */
  double cceom_energy;
  double R0;
  double L0;
  int L_irr;
  int R_irr;
  int G_irr;
  double L_root;
  double R_root;
  char *wfn;
};

