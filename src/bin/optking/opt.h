#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

#ifdef C_EXTERN
# undef C_EXTERN
# define C_EXTERN extern
#else
# define C_EXTERN
#endif

#define OPT_AUX (13)
#define SQR(A) ((A)*(A))
/*--- These are just internal type specifiers ---*/
#define STRE_TYPE (0)
#define BEND_TYPE (1)
#define TORS_TYPE (2)
#define OUT_TYPE (3)
#define LIN_BEND_TYPE (4)
#define NUM_INTCO_TYPES (4)
#define PRINT_TO_GEOM (113)
#define PRINT_TO_30 (114)
/* Limits to hard-wired arrays */
#define MAX_SALCS (500)
#define MAX_ATOMS (300)
#define MAX_ZVARS (500)
#define MAX_LINELENGTH (133)
#define MAX_SALC_LENGTH (1000)
#define MAX(I,J) ((I>J) ? I : J)
#define MIN(I,J) ((I>J) ? J : I)
#define EVAL_TOL (1.0E-14)                                     /* tolerance for eigenvalues (used in sq_rsp() and irrep() ) */
#define REDUNDANT_EVAL_TOL (1.0E-10)
#define SPANNED_IRREP_TOL (0.05)                               /* if character greater than this, irrep projected and kept */
#define LABEL_LENGTH (4) // for point group and irrep labels
/* step size limits */
#define STEP_LIMIT (0.1)     /* max step size if coord has small value */
#define STEP_PERCENT (0.3)   /* if coord large valued, max percentage allowed for step */
#define NONLINEAR_DIST (4.0E-3)
#define MIN_DQ_STEP (1.0E-12)
#define MIN_CART_OUT (1.0E-12)
#define MIN_LIN_COS (1.0E-10)

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN FILE *fp_input, *fp_intco, *fp_fconst, *fp_opt_aux, *fp_11;
EXTERN int *ops_in_class;
EXTERN int nirreps, *irr;
EXTERN int num_nonzero;      /* number of non-redundant di coordinates (eigenvectors of G with nonzero eigenvalues) */
EXTERN char ptgrp[4];        /*molecular point group*/
/* EXTERN int *number_internals; */
/* EXTERN double energy_chkpt(void); */

#ifdef C_CODE
C_EXTERN void punt(char *message);
C_EXTERN void open_PSIF(void);
C_EXTERN void close_PSIF(void);
C_EXTERN void exit_io(void);
C_EXTERN void print_mat2(double **matrix, int rows, int cols, FILE *of);
C_EXTERN void print_mat5(double **matrix, int rows, int cols, FILE *of);
C_EXTERN void cross_product(double *u,double *v,double *out);
C_EXTERN void scalar_mult(double a, double *vect, int dim);
C_EXTERN void scalar_div(double a, double *vect);
C_EXTERN int div_int(int big, int little);
C_EXTERN double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);
C_EXTERN double energy_chkpt(void);
C_EXTERN void dgeev_optking(int L, double **G, double *lambda, double **alpha);
C_EXTERN double **mass_mat(double *masses);
C_EXTERN double **unit_mat(int dim);
C_EXTERN void swap(int *a, int *b);
C_EXTERN void swap_tors(int *a, int *b, int *c, int *d);
C_EXTERN void zval_to_symbol(double zval, char *sym);

/* C_EXTERN int **get_char_table(char *ptgrp); returns the character table 
   C_EXTERN char **get_symm_ops(char *ptgrp); "     " symm operation labels 
   C_EXTERN int *get_ops_in_class(char *ptgrp);
   #include <dmalloc.h> */
#else
extern "C" void punt(char *message);
extern "C" void open_PSIF(void);
extern "C" void close_PSIF(void);
extern "C" void exit_io(void);
extern "C" void print_mat2(double **matrix, int rows, int cols, FILE *of);
extern "C" void print_mat5(double **matrix, int rows, int cols, FILE *of);
extern "C" void cross_product(double *u,double *v,double *out);
extern "C" void scalar_mult(double a, double *vect, int dim);
extern "C" void scalar_div(double a, double *vect);
extern "C" int div_int(int big, int little);
extern "C" double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);
extern "C" double energy_chkpt(void);
extern "C" void dgeev_optking(int L, double **G, double *lambda, double **alpha);
extern "C" double **mass_mat(double *masses);
extern "C" double **unit_mat(int dim);
extern "C" void swap(int *a, int *b);
extern "C" void swap_tors(int *a, int *b, int *c, int *d);
extern "C" void zval_to_symbol(double zval, char *sym);
/* extern "C" int **get_char_table(char *ptgrp); returns the character table 
   extern "C"  char **get_symm_ops(char *ptgrp); "     " symm operation labels
   extern "C" { #include <dmalloc.h> } */
#endif

#define MODE_DISP_NOSYMM   (10)
#define MODE_DISP_IRREP    (11)
#define MODE_DISP_LOAD     (12)
#define MODE_DISP_USER     (13)
#define MODE_LOAD_REF      (14)
#define MODE_OPT_STEP      (15)
#define MODE_FREQ_ENERGY   (16)
#define MODE_GRAD_ENERGY   (17)
#define MODE_FREQ_GRAD_NOSYMM (18)
#define MODE_FREQ_GRAD_IRREP  (19)
#define MODE_GRAD_SAVE        (20)
#define MODE_ENERGY_SAVE      (21)

struct OPTInfo {

  int mode;
  int disp_num;
  int points;
  int points_freq;
  int irrep;
  int simples_present;
  int salcs_present;
  int constraints_present;
  int nconstraints;
  int *constraints;

/* print options */
  int print_simples;
  int print_params;
  int print_delocalize;
  int print_symmetry;
  int print_hessian;

/* optimization parameters */
  int optimize;
  int redundant;
  int delocalize;
  int do_only_displacements;
  int zmat;
  int zmat_simples;
  int bfgs;
  int bfgs_use_last;
  int dertype;
  int numerical_dertype;
  int iteration;
  int micro_iteration;
  double conv; /* MAX force */
  double ev_tol;
  double scale_connectivity;
  double disp_size;
  int mix_types;
  int natom;
  int nallatom;
  int *atom_dummy;
  int *to_dummy;
  int *to_nodummy;
  int dummy_axis_1;
  int dummy_axis_2;

/* Back-transformation parameters */
  int bt_max_iter;
  double bt_dq_conv;
  double bt_dx_conv;

/* Obscure limits in intco evaluation */
  double cos_tors_near_1_tol;
  double cos_tors_near_neg1_tol;
  double sin_phi_denominator_tol;
};

struct SYMInfo {
  char *symmetry;
  int nirreps;
  int **ct;
  char **irrep_lbls;
  char **clean_irrep_lbls;
  char **op_lbls;
  int **ict;
  int **fict;
  int **ict_ops;
  int **ict_ops_sign;
};

EXTERN struct OPTInfo optinfo;
EXTERN struct SYMInfo syminfo;









