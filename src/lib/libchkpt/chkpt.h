#include <chkpt_params.h>

#define MAX_ELEMNAME 13

/* Uncomment after removal of old file30 */
/*structure to hold z-matrix info*/
/*  struct z_entry { */
/*    int bond_atom;   */         /*first reference atom (bond)*/
/*    int angle_atom;    */       /*second reference atom (angle)*/
/*    int tors_atom;    */        /*third reference atom (torsion)*/
/*    int bond_opt;    */         /*flags indicating to optimize values*/
/*    int angle_opt; */
/*    int tors_opt; */
/*    double bond_val; */          /*coordinate values*/
/*    double angle_val; */
/*    double tors_val; */
/*    char bond_label[20]; */     /*variable labels, if any*/
/*    char angle_label[20]; */
/*    char tors_label[20]; */
/*    }; */

/*--- Types of reference determinants ---*/
/* Uncomment after removal of old file30 */
/* typedef enum {ref_rhf = 0, ref_uhf = 1, ref_rohf = 2, ref_tcscf = 3,
   ref_rks = 4, ref_uks = 5} reftype; */


int chkpt_init(void);
int chkpt_close(void);

int chkpt_rd_ncalcs(void);

double chkpt_rd_escf(void);
void chkpt_wt_escf(double);

double chkpt_rd_etot(void);
void chkpt_wt_etot(double);

int chkpt_rd_phase_check(void);
void chkpt_wt_phase_check(int);

int chkpt_rd_iopen(void);
void chkpt_wt_iopen(int);

int chkpt_rd_ref(void);
void chkpt_wt_ref(int);

int chkpt_rd_nmo(void);
void chkpt_wt_nmo(int);

int chkpt_rd_nsymhf(void);
void chkpt_wt_nsymhf(int);

int chkpt_rd_mxcoef(void);
void chkpt_wt_mxcoef(int);

int *chkpt_rd_atom_position(void);
void chkpt_wt_atom_position(int *);

int **chkpt_rd_ict(void);
void chkpt_wt_ict(int **);

double *chkpt_rd_zvals(void);
void chkpt_wt_zvals(double *zvals);

double *chkpt_rd_exps(void);
void chpt_wt_exps(double *);

int *chkpt_rd_us2s(void);
void chkpt_wt_us2s(int *);

int *chkpt_rd_orbspi(void);
void chkpt_wt_orbspi(int *);

int *chkpt_rd_clsdpi(void);
void chkpt_wt_clsdpi(int *);

int *chkpt_rd_openpi(void);
void chkpt_wt_openpi(int *);

int *chkpt_rd_sopi(void);
void chkpt_wt_sopi(int *);

char *chkpt_rd_label(void);
void chkpt_wt_label(char *);

char **chkpt_rd_irr_labs(void);
void chkpt_wt_irr_labs(char **);

double *chkpt_rd_contr(void);
void chkpt_wt_contr(double *);
double **chkpt_rd_contr_full(void);

int *chkpt_rd_sprim(void);
void chkpt_wt_sprim(int *);

int *chkpt_rd_snuc(void);
void chkpt_wt_snuc(int *);

int *chkpt_rd_stype(void);
void chkpt_wt_stype(int *);

int *chkpt_rd_snumg(void);
void chkpt_wt_snumg(int *);

int *chkpt_rd_sloc(void);
void chkpt_wt_sloc(int *);

int **chkpt_rd_shell_transm(void);
void chkpt_wt_shell_transm(int **);

int chkpt_rd_nentry(void);
void chkpt_wt_nentry(int);

char **chkpt_rd_fement(void);
void chkpt_wt_felement(char **);

double **chkpt_rd_usotao(void);
void chkpt_wt_usotao(double **);

double **chkpt_rd_usotbf(void);
void chkpt_wt_usotbf(double **);

int *chkpt_rd_sloc_new(void);
void chkpt_wt_sloc_new(int *);

int *chkpt_rd_ua2a(void);
void chkpt_wt_ua2a(int *);

int *chkpt_rd_symoper(void);
void chkpt_wt_symoper(int *);

char *chkpt_rd_sym_label(void);
void chkpt_wt_sym_label(char *sym_label);

struct z_entry *chkpt_rd_zmat(void);
void chkpt_wt_zmat(struct z_entry *);

int *chkpt_rd_shells_per_am(void);
void chkpt_wt_shells_per_am(int *);

int *chkpt_rd_am2canon_shell_order(void);
void chkpt_wt_am2canon_shell_order(int *);

double **chkpt_rd_rref(void);
void chkpt_wt_rref(double **);

double **chkpt_rd_fgeom(void);
void chkpt_wt_fgeom(double **);

double **chkpt_rd_geom(void);
void chkpt_wt_geom(double **);

double chkpt_rd_enuc(void);
void chkpt_wt_enuc(double);

int chkpt_rd_num_unique_atom(void);
void chkpt_wt_num_unique_atom(int);

int chkpt_rd_num_unique_shell(void);
void chkpt_wt_num_unique_shell(int);

int chkpt_rd_rottype(void);
void chkpt_wt_rottype(int);

int chkpt_rd_max_am(void);
void chkpt_wt_max_am(int);

int chkpt_rd_nso(void);
void chkpt_wt_nso(int);

int chkpt_rd_nao(void);
void chkpt_wt_nao(int);

int chkpt_rd_nshell(void);
void chkpt_wt_nshell(int);

int chkpt_rd_nirreps(void);
void chkpt_wt_nirreps(int);

int chkpt_rd_nprim(void);
void chkpt_wt_nprim(int);

int chkpt_rd_natom(void);
void chkpt_wt_natom(int);

double *chkpt_rd_evals(void);
double *chkpt_rd_alpha_evals(void);
double *chkpt_rd_beta_evals(void);
void chkpt_wt_evals(double *);
void chkpt_wt_alpha_evals(double *);
void chkpt_wt_beta_evals(double *);

double **chkpt_rd_scf(void);
double **chkpt_rd_alpha_scf(void);
double **chkpt_rd_beta_scf(void);
void chkpt_wt_scf(double **);
void chkpt_wt_alpha_scf(double **);
void chkpt_wt_beta_scf(double **);

double **chkpt_rd_scf_irrep(int);
double **chkpt_rd_alpha_scf_irrep(int);
double **chkpt_rd_beta_scf_irrep(int);
void chkpt_wt_scf_irrep(double **, int);
void chkpt_wt_alpha_scf_irrep(double **, int);
void chkpt_wt_beta_scf_irrep(double **, int);

double **chkpt_rd_lagr(void);
double **chkpt_rd_alpha_lagr(void);
double **chkpt_rd_beta_lagr(void);
void chkpt_wt_lagr(double **);
void chkpt_wt_alpha_lagr(double **);
void chkpt_wt_beta_lagr(double **);

double **chkpt_rd_ccvecs(void);
void chkpt_wt_ccvecs(double **);

double **chkpt_rd_schwartz(void);
void chkpt_wt_schwartz(double **);

double chkpt_rd_ecorr(void);
void chkpt_wt_ecorr(double);

double chkpt_rd_eref(void);
void chkpt_wt_eref(double);

double chkpt_rd_efzc(void);
void chkpt_wt_efzc(double);

int chkpt_rd_disp(void);
void chkpt_wt_disp(int);
