/* Input parameters for ccenergy */
struct Params {
  long int memory;
  int aobasis;
  int cachelev;
  int cachetype;
  int ref;
  int eom_ref;
  int local;
};

struct Eom_params {
  int max_iter;
  int vectors_per_root;
  int *states_per_irrep;
  int *cs_per_irrep;
  double eval_tol;
  double residual_tol;
  int prop_root;
  int prop_sym;
  int save_all;
  int print_singles;
  double complex_tol;
  double schmidt_add_residual_tol;
  int max_iter_SS;
  int vectors_per_root_SS;
  int excitation_range;
  double residual_tol_SS;
  char *guess;
  int dot_with_Lg; // check orthogonality with Lg (must run cclambda first)
};

