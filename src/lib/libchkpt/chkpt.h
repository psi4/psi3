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

int chkpt_init(void);
int chkpt_close(void);

int chkpt_rd_iopen(void);
int chkpt_rd_max_am(void);
int chkpt_rd_mxcoef(void);
int chkpt_rd_nao(void);
int chkpt_rd_natom(void);
int chkpt_rd_ncalcs(void);
int chkpt_rd_nirreps(void);
int chkpt_rd_nmo(void);
int chkpt_rd_nprim(void);
int chkpt_rd_nshell(void);
int chkpt_rd_nso(void);
int chkpt_rd_nsymhf(void);
int chkpt_rd_num_unique_atom(void);
int chkpt_rd_num_unique_shell(void);
int chkpt_rd_phase_check(void);
int chkpt_rd_ref(void);
int chkpt_rd_rottype(void);
int *chkpt_rd_am2canon_shell_order(void);
int *chkpt_rd_atom_position(void);
int *chkpt_rd_shells_per_am();
int *chkpt_rd_sloc();
int *chkpt_rd_sloc_new();
int *chkpt_rd_snuc();
int *chkpt_rd_snumg();
int *chkpt_rd_sprim();
int *chkpt_rd_sopi();
int *chkpt_rd_stype();
int *chkpt_rd_symoper();
int *chkpt_rd_ua2a();
int *chkpt_rd_us2s();
int *chkpt_rd_orbspi(void);
int *chkpt_rd_clsdpi(void);
int *chkpt_rd_openpi(void);
int *chkpt_rd_sopi(void);

char *chkpt_rd_label(void);
char **chkpt_rd_irr_labs(void);
char **chkpt_rd_hfsym_labs(void);
