/**************************
  Prototypes for functions
 **************************/

/*
   strange, but necessary, gprgid function
*/
char *gprgid ();


/*
   Function to perform symmetry operations on input cartesians
*/
int operation(int op_num, char *A, double **Rotate, double *B, double *D);

/* 
   C2z function
*/
void c2z(double *A, double *B);

/* 
   C2y function 
*/
void c2y(double *A, double *B);

/* 
   C2x function 
*/
void c2x(double *A, double *B);

/* 
   inversion function 
*/
void inversion(double *A, double *B);

/* 
   sig_xy function 
*/
void sig_xy(double *A, double *B);

/* 
   sig_xz function 
*/
void sig_xz(double *A, double *B);

/* 
   sig_yz function 
*/
void sig_yz(double *A, double *B);

/* 
   quit routine 
*/
void punt(const char *msg);

/*
   Function to get Atomic Number
*/
void atom_num(char *A,  double *C);

/*
   Function to calculate the unit vectors between each pair of atoms and
   store them in E. 
*/
void unit_vectors(int num_atoms, double *X, double *Y, double *Z, 
                  double **Distance, double ***E);

/*
   Function to compute the unit vector between a pair of atoms
*/
void unit_vec(double *, double *, double *);

/*
   Function to calculate all non-redundant bond angles
*/
void calc_bond_angles(double ***E, 
                      double ***Bond_Angle, 
                      double **Distance, FILE *outfile);

/*
   Takes a dot product
*/
double dot_prod(double *A, double *B);

/*
   Computes the cross product
*/
void cross_prod(double *v1, double *v2, double *out);


/*
   Function to calculate all internuclear distances
*/
void calc_distance(double **, double *, int);

/*
   Function to calculate the nuclear repulsion energy
*/
void Nuc_repulsion(double *Distance, double *repulsion);

/* 
   Create a rotation matrix - used when dropping symmetry
*/
void Rotate(double **Rotate_axis, double *old_coord);

/*
  Function rotates geometry[][] by new_coord[][]
*/

void rotate_geometry(double **geom, double **new_coord);

/*
  Function rotates full_geom[][] by new_coord[][]
*/

void rotate_full_geom(double **geom, double **new_coord);

/*
  Function updates Rref to remember the effect of the rotation
  described by R
*/

void memorize_rotation(double **R);

/*
   Convert from angstroms to bohr
*/
void conv_2_bohr(double *old_coord);

/*
   Create the elem_name matrix
*/
void init_elem_names();

/*
   Main routine to read the basis set info.  It calls recur, and all of the normalization
   routines. */
void read_basis();

/*
   normalization functions in read_basis.c 
*/
void setup();
void normalize();
int parse_am();
double int_pow();
double ovlp();


/*
   This function calls itself recursively while reading the basis set info. in 
   pbasis.dat (or in input deck) until it has gotten past all of the "GET" levels */
void recur(char **ip_token1, char **ip_token2, int num_levels, 
int atom_number, double **basis_set,
int *count1, int *count2, int num_exponents, int *fprim, int* lprim,
int *ang_mom);


char *init_char_array(int B);
char **init_char_matrix(int A, int B);

void start_io();
void stop_io();

void init_oldcalc();
void cleanup_oldcalc();
void oldcalc_projection();

void init_gto(int);
void cleanup_gto(int);
void save_oldmos();

double **overlap_new_old();
double **overlap();
void write_scf_to_file30();

void parsing();
void parsing_cmdline(int argc, char *argv[]);

/*
  Make canonical and reference frames equivalent
  */
void canon_eq_ref_frame();

