/***************************
    Global variables
 ***************************/

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

/*--- NOTE!
  Whenever I refer to basis function - it means either puream/cart function
  AO means cartesian Gaussian
  ---*/


/*Super-global stuff - the same no matter what calculation is running */
EXTERN FILE *infile, *outfile;
EXTERN int *ioff;
EXTERN int *df;                     /*df[i] = (i-1)!!*/
EXTERN char **elem_name;            /*Element names*/

/*Calculation options*/
EXTERN int cartOn;                  /*Cartesian input flag*/
EXTERN int puream;
EXTERN int shownorm;		    /*Show normalized basis set*/
EXTERN char *units;
EXTERN int print_lvl;               /*Printing level*/
EXTERN int no_reorient;             /*No reorientation?*/


/*Labels*/
EXTERN char *label;                 /*Label*/


/*Symmetry-related stuff*/
EXTERN char *symmetry;              /*Symmetry group label*/
EXTERN char *subgroup;              /*Subgroup label*/
EXTERN char *unique_axis;           /*Unique axis label*/
EXTERN int **irr_char;              /*Character table for irreducible representations of D2h and subgroups*/
EXTERN int *irr_char_str;           /*Character table for irreps packed in bytes: "-1" = 1, "1" = 0*/
EXTERN char **irr_labels;           /*Labels for irreducible representations*/
EXTERN int *sym_oper;               /*Array that maps symmetry operation number of the point group to the canonical
				      operation numbering (see symmetry.h)*/


/*Calculation-dependent scalars*/
EXTERN double conv_factor;          /*Conversion factor for geometry:
				      if UNITS=BOHR or AU - 1.0,
				      if UNITS=ANGSTROM or ANGSTROMS - 1/_bohr2angstroms.*/
typedef enum {                      /*Type of the rigid rotor (3 + number of zero moments - number of */
                                    /*non-zero unique moments) */
    asymmtop = 0,                   /*0 - asymmetric top (3 + 0 - 3) */
    symmtop = 1,                    /*1 - symmetric top (3 + 0 - 2) */
    sphtop = 2,                     /*2 - spherical top (3 + 0 - 1) */
    linear = 3,                     /*3 - linear molecule (3 + 1 - 1)*/
    atom = 6                        /*6 - atom (3 + 3 - 0) */
} rotortype;
EXTERN rotortype rotor;             /*Type fo the rotor we are dealing with*/
EXTERN int nirreps;                 /*Number of irreducible representations in the computational (largest Abelian)
				      point group*/
EXTERN int num_atoms;               /*Total number of atoms*/
EXTERN int num_uniques;             /*Number of unique atoms*/
EXTERN int num_shells;              /*Total number of shells*/
EXTERN int num_unique_shells;       /*Number of unique shells*/
EXTERN int num_ao;                  /*Total number of AOs*/
EXTERN int num_so;                  /*Total number of SOs*/
EXTERN int num_prims;               /*Number of unique primitives*/
EXTERN int max_angmom;              /*Maximum angular momentum type in the basis*/
EXTERN int num_classes;             /*Number of atom classes*/
EXTERN int num_unique_classes;      /*Number of symmetry unique atom classes*/
EXTERN int num_angso_coeff;
EXTERN int ap_flag;         /*Bitfield of flags indicating presence of atoms in certain positions
			      bit 1 = on C2x axis (position = 2)
			      bit 2 = on C2y axis (position = 4)
			      bit 3 = on C2z axis (position = 8)
			      bit 4 = in inversion center (position = 16)
			      bit 5 = in sig_xy plane (position = 32)
			      bit 6 = in sig_xz plane (position = 64)
			      bit 7 = in sig_yz plane (position = 128)
			      bit 0 = in general position (position = 1) */

/*Calculation-dependent arrays*/
EXTERN double **geometry;	    /*Cartesian coordinates (in a.u.) of atoms*/
EXTERN double *nuclear_charges;	    /*Nuclear charges*/
EXTERN char **element;       	    /*Atom names*/
EXTERN char **atom_basis;           /*Array of basis set names*/
EXTERN int **atom_orbit;            /*Atom orbits*/
EXTERN int **class_orbit;           /*Class orbits*/
EXTERN int **red_unique_orbit;      /*Reduced atomic orbits for unique atoms (each symmetry equiv atom appears only once),
				      unique atom itself is included*/
EXTERN int *u2a;                    /*Mapping of unique atom number on the full atom list*/
EXTERN int *uc2c;                   /*Mapping of unique class number on the full class list*/
EXTERN int *unique_degen;           /*Degeneracy of unique atoms*/
EXTERN int *unique_class_degen;     /*Degeneracy of unique classes*/


/*Irrep-dependent arrays*/
EXTERN int *num_cart_so_per_irrep;  /*Number of cartesian SOs per irrep*/
EXTERN int *num_so_per_irrep;       /*Number of cart/pureang SOs per irrep*/
EXTERN int **num_cart_so;           /*Number of cartesian type SOs in each irrep for each angular momentum type*/
EXTERN int **num_pureang_so;        /*Number of pure ang. momentum type SOs in each irrep for each angular momentum type*/
EXTERN int **num_redun_so;          /*Difference between num_cart_so and num_pureang_so*/


/*Basis set arrays - written to file30*/
EXTERN int *nshells_per_atom;         /*Number of shells per atom*/
EXTERN int *first_shell_on_atom;      /*Number of the first shell from an atom*/
EXTERN double *exponents;             /*Exponents*/
EXTERN double *contr_coeff;           /*Contraction coefficients*/
EXTERN int *first_prim_shell;         /*Number of the first primitive for a shell*/
EXTERN int *shell_nucleus;            /*Number of the nucleus a shell belongs to*/
EXTERN int *shell_ang_mom;            /*Angular momentum of a shell*/
EXTERN int *nprim_in_shell;           /*Number of primitives in a shell*/
EXTERN int *first_ao_shell;           /*Number of the first AO from a shell*/
EXTERN int *first_basisfn_shell;      /*Number of the first basis function from a shell*/
EXTERN int *first_ao_type_shell;      /*Type(xx,yy,etc.) of the first AO from a shell*/
EXTERN int *last_ao_type_shell;       /*Type(xx,yy,etc.) of the last AO from a shell */
EXTERN double *angso_coeff;           /*Coefficients of SO's*/
EXTERN int **angso_labels;            /*Labels of SO matrices*/
EXTERN int **first_angso_coeff_class; /*Pointer to the first coefficient in angso_coeff for nth class and angmom l*/
EXTERN int **last_angso_coeff_class;  /*Pointer to the last coefficient in angso_coeff for nth class and angmom l*/
EXTERN double ***ao_type_transmat;    /*Transformation matrices for AO types*/
EXTERN int *pureang_so_l;             /*Angular momentum number of a pure angular momentum type SO*/
EXTERN int *pureang_so_m;             /*Modulus of m of a pure angular momentum type SO*/


/*SO to AO transformation-related stuff*/
EXTERN double ***cart2pureang;             /*Basic cartesian to pure angular momentum matrices*/
EXTERN double **usotao_cart;               /*AO to cartesian SO matrix*/
EXTERN double **cart_usotao_pureang;       /*Cartesian SO to pure angular momentum SO matrix*/
EXTERN double **usotao;                    /*AO to SO matrix*/
EXTERN double **usotbf;                    /*Basis functions to SO matrix*/


/*Symmetry-class-related information*/
EXTERN int *max_angmom_class;         /*Maximum angular momentum of basis functions on atoms of nth class*/
EXTERN int *atom_class;               /*Class an atom belongs to*/
EXTERN int *atom_position;            /*Symmetry position of an atom*/
EXTERN double ****class_so_coeff;          /*Matrix of AO to cartesian SO for each class and angular momentum type*/
EXTERN int ***num_cart_so_in_class;        /*Number of cartesian SO for each class, ang. momentum type and irrep*/
EXTERN int ***num_pureang_so_in_class;     /*Number of pure angular momentum SO for each class, ang. momentum type and irrep*/


/*Auxiliary arrays*/
EXTERN int **ao_type_irr;        /*Irreducible representation an AO of a given type positioned in the origin belongs to*/


/* Arrays of x, y, and z exponents in cartesian Gaussians (AOs)*/
EXTERN int **xexp_ao, **yexp_ao, **zexp_ao;

