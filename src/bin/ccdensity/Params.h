/* Input parameters for cclambda */
struct Params {
  double tolerance;
  long int memory;
  int cachelev;
  int aobasis;
  int ref;
  int onepdm; /* produce ONLY the onepdm for properties */
  int relax_opdm;
  int use_zeta;
  int calc_xi;
  int connect_xi;
  int restart;
  int ground;
  int transition; 
  int dertype;
  double cceom_energy;
  double R0;
  double L0;
  int L_irr;
  int R_irr;
  int G_irr;
  int L_root;
  int R_root;
  char *wfn;
  double overlap1; /* <L1|R1> */
  double overlap2; /* <L2|R2> */
  double RD_overlap; /* Rmnef <mn||ef> */
  double RZ_overlap; /* <R|zeta> */
  char *gauge;
  int nstates;
  int prop_sym;
  int prop_root;
};

struct TD_Params {
  int irrep;
  int root;
  double R0;
  double cceom_energy;
  char L1A_lbl[32];
  char L1B_lbl[32];
  char L2AA_lbl[32];
  char L2BB_lbl[32];
  char L2AB_lbl[32];
  char R1A_lbl[32];
  char R1B_lbl[32];
  char R2AA_lbl[32];
  char R2BB_lbl[32];
  char R2AB_lbl[32];
  double OS;
  double RS;
};
