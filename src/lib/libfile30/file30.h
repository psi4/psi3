#ifndef FILE30_H
#define FILE30_H

#include <file30_params.h>

/*structure to hold z-matrix info*/
struct z_entry {
  int bond_atom;            /*first reference atom (bond)*/
  int angle_atom;           /*second reference atom (angle)*/
  int tors_atom;            /*third reference atom (torsion)*/
  int bond_opt;             /*flags indicating to optimize values*/
  int angle_opt;
  int tors_opt;
  double bond_val;          /*coordinate values*/
  double angle_val;
  double tors_val;
  char bond_label[20];      /*variable labels, if any*/
  char angle_label[20];
  char tors_label[20];
  };

/*--- Types of reference determinants ---*/
typedef enum {ref_rhf = 0, ref_uhf = 1, ref_rohf = 2, ref_tcscf = 3,
      ref_rks = 4, ref_uks = 5} reftype;

int file30_init(void);
int file30_close(void);
double **file30_rd_rref(void);
double **file30_rd_schwartz(void);
double **file30_rd_alpha_lagr(void);
double **file30_rd_beta_lagr(void);
double **file30_rd_lagr(void);
double **file30_rd_usotbf(void);
double **file30_rd_usotao_new(void);
double **file30_rd_contr_full(void);
double **file30_rd_alpha_blk_scf(int);
double **file30_rd_beta_blk_scf(int);
double **file30_rd_blk_scf(int);
double **file30_rd_ccvecs(void);
double **file30_z_to_cart(struct z_entry *z_geom, int num_atoms);
double *file30_rd_contr(void);
double file30_rd_eref(void);
double file30_rd_ecorr(void);
double file30_rd_efzc(void);
double file30_rd_enuc(void);
double file30_rd_escf(void);
double **file30_rd_fgeom(void);
double *file30_rd_alpha_evals(void);
double *file30_rd_beta_evals(void);
double *file30_rd_evals(void);
double *file30_rd_exps(void);
double **file30_rd_geom(void);
double **file30_rd_alpha_scf(void);
double **file30_rd_beta_scf(void);
double **file30_rd_scf(void);
struct z_entry *file30_rd_zmat(void);
double *file30_rd_zvals(void);
char *file30_rd_sym_label(void);
char *file30_rd_corr_lab(void);
char **file30_rd_felement(void);
char **file30_rd_hfsym_labs(void);
char **file30_rd_irr_labs(void);
char *file30_rd_label(void);
int *file30_rd_shells_per_am(void);
int *file30_rd_am2canon_shell_order(void);
int *file30_rd_ua2a(void);
int *file30_rd_symoper(void);
int *file30_rd_sloc_new(void);
int *file30_rd_us2s(void);
int *file30_rd_atom_position(void);
int *file30_rd_sopi(void);
int *file30_rd_clsdpi(void);
int **file30_rd_ict(void);
int **file30_rd_shell_transm(void);
int file30_rd_max_am(void);
int file30_rd_rottype(void);
int file30_rd_num_unique_atom(void);
int file30_rd_num_unique_shell(void);
int file30_rd_iopen(void);
int file30_rd_ref(void);
int file30_rd_mxcoef(void);
int file30_rd_nao(void);
int file30_rd_natom(void);
int file30_rd_ncalcs(void);
int file30_rd_nentry(void);
int file30_rd_nirreps(void);
int file30_rd_nmo(void);
int file30_rd_nso(void);
int file30_rd_nprim(void);
int file30_rd_nshell(void);
int file30_rd_nsymhf(void);
int file30_rd_ref(void);
int *file30_rd_openpi(void);
int *file30_rd_orbspi(void);
int file30_rd_phase_check(void);
int *file30_rd_scf_ptrs(void);
int *file30_rd_sloc(void);
int *file30_rd_smax(void);
int *file30_rd_smin(void);
int *file30_rd_snuc(void);
int *file30_rd_snumg(void);
int *file30_rd_sprim();
int *file30_rd_stype();
int file30_rd_disp();
void pack_4int(int **, int *, int, int);
void unpack_4int(int *, int **, int, int);
void file30_wt_eref(double);
void file30_wt_efzc(double);
void file30_wt_enuc(double);
void file30_wt_shell_transm(int **);
void file30_wt_isc(int **);
void file30_wt_ipc(int **);
void file30_wt_ecorr(double);
void file30_wt_escf(double);
void file30_wt_fgeom(double **);
void file30_wt_clsdpi(int *);
void file30_wt_corr_lab(char *);
void file30_wt_geom(double **);
void file30_wt_iopen(int);
void file30_wt_openpi(int *);
void file30_wt_alpha_scf(double **);
void file30_wt_beta_scf(double **);
void file30_wt_scf(double **);
void file30_wt_alpha_evals(double *);
void file30_wt_beta_evals(double *);
void file30_wt_evals(double *);
void file30_wt_zvals(double *);
void file30_wt_alpha_blk_scf(double **,int);
void file30_wt_beta_blk_scf(double **,int);
void file30_wt_blk_scf(double **,int);
void file30_wt_zmat(struct z_entry *z_geom, int num_atoms);
void file30_wt_disp(int disp);

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#endif
