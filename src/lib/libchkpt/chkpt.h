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

int chkpt_rd_iopen();
int chkpt_rd_max_am();
int chkpt_rd_mxcoef();
int chkpt_rd_nao();
int chkpt_rd_natom();
int chkpt_rd_ncalcs();
int chkpt_rd_nirreps();
int chkpt_rd_nmo();
int chkpt_rd_nprim();
int chkpt_rd_nshell();
int chkpt_rd_nso();
int chkpt_rd_nsymhf();
int chkpt_rd_num_unique_atom();
int chkpt_rd_num_unique_shell();
int chkpt_rd_phase_check();
int chkpt_rd_ref();
int chkpt_rd_rottype();
int *chkpt_rd_am2canon_shell_order();
int *chkpt_rd_atom_position();
