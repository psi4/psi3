/*-----------------------
  structure declarations
 -----------------------*/

/*- 2-index quantity -*/
struct oebuf {
  double val;
  int ij;
};

/*- 4-index quantity with index range of (0,65535) -*/
struct tebuf {
  double val;
  short int i;
  short int j;
  short int k;
  short int l;
};


struct gaussian_function{
   double exp;		/* orbital exponent */
   double ccoeff[MAX_AM];	/* comb. of contraction coeff and normalization */
};

struct shell_def{
   int center;		/* atom on which shell is centered */
   int am;		/* angular momentum of shell (+1?)*/
   int n_prims;		/* number of primitives in shell */
   int fprim;		/* pointer to the first primitive in shell */
   int fbf;             /* pointer to the first basis function from shell */
   int fao;             /* pointer to the first AO from shell */
   int *trans_vec;      /* shell symmetry transformation vector */
};

struct coordinates{
   double x;  /*  what do you think these are? */
   double y;
   double z;
   double Z_nuc; /* nuclear charge */
};

struct shell_pair{
  int i, j;
  double ***P;
  double AB[3];
  double ***PA;
  double ***PB;
  double *a1, *a2, **gamma;
  double *inorm, *jnorm;
  double **Sovlp;
  /*--- Having open-shell and closed-shell and alpha and beta densities
        in UHF case is redundant. Thus only one pair is used at a given time, and the
	other pair is empty (ptr=NULL)!

	In fact, the difference open- and closed-shell densities in UHF case are
	used in hf_fock to form the Fock matrix incrementally.
	the total alpha and beta densities in UHF case are used in DFT to form the total
	XC matrix. ---*/
  double **dmat;
  double **dmato;
  double **dmatb;
  double **dmata;
  /*--- SCF Lagrangian ---*/
  double **lagr;
  double Smax;
  double Dmax;
};

struct unique_shell_pair{
  int *SOpair_npi;
  int **SOpair_so_i;
  int **SOpair_so_j;
  int **SOpair_bf_i;
  int **SOpair_bf_j;
};

/*--- Matrix element (to be used to handle sparse matrices) ---*/
typedef struct {
  int row;
  int column;
  double value;
} mat_elem;


/*--- These are basically identical, left from Justin ---*/
struct double_array {
 int n1;
 double *d;
};
typedef struct double_array double_array_t;

struct struct_double_matrix {
  int n1;
  int n2;
  double **d /*[n1]*/ /*[n2]*/;
  };
typedef struct struct_double_matrix double_matrix_t;


enum scftype {rhf = 0, uhf = 1, rohf = 2, twocon = 3};

typedef struct {
    char *wfn;                         /* Wavefunction */
    char *dertype;                     /* Derivative type */
    double cutoff;                     /* Cutoff on ERIs/Fock matrix elements */
    double hf_exch;                    /* Portion of exact HF exchange in the Fock matrix */
    int print_lvl;                     /* Print level */
    int max_memory;                    /* Maximum amount of memory to use, in double words */
    int memory;                        /* Amount left available */
    int make_oei;                      /* Flag to compute one-electron integrals */
    int make_fock;                     /* Flag to compute Fock matrix */
    int make_eri;                      /* Flag to compute two-electron integrals */
    int make_deriv1;                   /* Flag to compute first derivatives of one- and two-electron integrals */
    int make_oeprop;                   /* Flag tp compute one-electron property integrals */
    int make_mp2;                      /* Flag to compute MP2 energy directly */
    int make_r12ints;                  /* Compute integrals for linear R12 methods */
    int make_mp2r12;                   /* Flag to compute MP2-R12 ebergy directly */
    int symm_ints;                     /* This flag should be set whe individual integrals over SO need to be computed */
    int scf_only;                      /* Means that ERIs will be used only in SCF calculations
					  (may save some space) */
    int num_threads;                   /* Number of threads */
    enum scftype reftype;              /* Reference type, e.g. RHF, ROHF, UHF */
} UserOptions_t;

typedef struct {
    int num_prims;                     /* number of primitive gaussians */
    int num_shells;                    /* number of shells */
    int max_num_prims;                 /* maximum number of primitives per shell */
    int puream;                        /* pure angular momentum flag */
    int num_ao;                        /* number of AO's */
    int max_am;                        /* maximum angular momentum in the basis */
    double **schwartz_eri;             /* the matrix num_shells by num_shells:
					  [si][sj] = max(ij|ij) i in si, j in sj  */
    struct shell_def *shells;          /* shell info */
    struct gaussian_function *cgtos;   /* cartesian gaussian information */
    struct shell_pair **shell_pairs;   /* shell pair info */
} BasisSet_t;

typedef struct {
    int num_unique_shells;             /* number of symmetry unique shells */
    int num_so;                        /* number of SO's */
    int nirreps;                       /* number of irreps */
    int max_stab_index;                /* maximum stabilizer index */
    int *atom_positions;               /* symmetry positions/stabilizers of atoms */
    int *us2s;                         /* unique shell number to full shell number mapping array */
    int *sopi;                         /* number of SO per irrep */
    int *sym_oper;                     /* mapping array between "canonical" and symmetry.h-defined
					  ordering of symmetry operations */
    int *so2symblk;                    /* SO number to symmetry block mapping array */
    int **dp_table;                    /* direct product multiplication table */
    int ***dcr;                        /* double coset representatives */
    int **dcr_dim;                     /* dimensions of double coset representatives */
    int **dcr_deg;
    int **GnG;
    char *symlabel;                    /* symmetry label */
    char **irr_labels;                 /* labels of irreps */
    double **usotao;                   /* SO to (basis functions if puream && !make_fock, AO otherwise)
					  transformation matrix */
    struct unique_shell_pair **us_pairs; /* unique shell symmetry info */
} SymmetryInfo_t;


typedef struct {
    int num_atoms;                     /* number of atoms */
    double Enuc;                       /* nuclear repulsion energy */
    char *label;                       /* calculation label */
    struct coordinates *centers;       /* nuclear centers */
} Molecule_t;
    
typedef struct {
    int itap30;               /* Checkpoint file */
    int itap33;               /* SO ERI file in IWL format */
    int itapS;                /* SO Overlap IWL file */
    int itapT;                /* SO Kinetic energy IWL file */
    int itapV;                /* SO Potential energy IWL file */
    int itapS_AO;             /* AO Overlap IWL file */
    int itapMX_AO;            /* AO mu(x) IWL file */
    int itapMY_AO;            /* AO mu(y) IWL file */
    int itapMZ_AO;            /* AO mu(z) IWL file */
    int itapDSCF;             /* "Interface" file between DSCF and CINTS */
    int itapD;                /* Correlated AO OPDM and Lagrangian from transqt */
    int itapG;                /* Correlated AO TPDM from transqt */
    int itapR12;              /* SO integrals of r12 operator */
    int itapT1;               /* SO integrals of [r12,T1] operator */
    int itapERI_MO;           /* MO ERI integrals */
    int itapR12_MO;           /* MO R12 integrals */
    int itapR12T2_MO;         /* MO [r12,T2] integrals */
} IOUnits_t;

typedef struct {
    double **bf_norm;                  /* normalization constants for cartesian GTOs of each
					 angular momentum level */
    double ***cart2pureang;            /* cartesian to pure angular momentum transformation matrices */
    double ****cc2pp;                  /* composite (CxC) cartesian to pure angular momentum transformation
					  matrices */
    mat_elem ****cc2pp_sparse;         /* sparse representation (row-compressed) of cc2pp */
    mat_elem ****pp2cc_sparse;         /* sparse representation (row-compressed) of the reverse of cc2pp */
    int ***cc2pp_rowlength;            /* this holds lengths of "compressed" rows in cc2pp_sparse */
    int ***pp2cc_rowlength;            /* see above */
} GTOs_t;

typedef struct {
    double Escf;              /* SCF energy */
    double Ecorr;             /* Correlation energy */
    double Eref;              /* Reference energy (if not SCF reference) */
    double *scf_evals[2];     /* SCF eigenvalues (alpha and beta spin) */
    double *scf_evals_occ[2]; /* Eigenvalues for active occupied orbitals in QTS order */
    double *scf_evals_uocc[2];/* Eigenvalues for active virtual orbitals in QTS order */
    double **scf_evec[2];     /* SCF eigenvectors in AO basis (alpha and beta spin)
			         NOTE: MOs are arranged in rows!!!!! */
    double **scf_evec_occ[2]; /* SCF eigenvectors in AO basis for all
				 doubly-occupied MOs in QTS order:
				 frozen DOCC MOs for each symmetry block come first,
				 then active DOCC MOs for each symmetry block */
    double **scf_evec_uocc[2];/* SCF eigenvectors in AO basis for all
				 vacant MOs in QTS order:
				 active UOCC MOs for each symmetry block come first,
				 then frozen UOCC MOs for each symmetry block */
    double tcscf_occ[2];      /* Squared coefficients of determinants in TCSCF wavefunction */
    double **Alpha, **Beta;   /* Alpha and Beta energy coupling coeffcients */
    int *mo2symblk;           /* Array that maps MO Pitzer index to its symblk number;
				 useful in manipulating Pitzer-indexed MOs */
    int *mo2symblk_occ[2];    /* Array that maps docc index to its symblk number;
			         useful in manipulating QTS-indexed MOs */
    int *mo2symblk_uocc[2];   /* Array that maps uocc index to its symblk number;
			         useful in manipulating QTS-indexed MOs */
    int *orbspi;              /* number of MOs per irrep */
    int *clsdpi;              /* number of doubly-occupied MOs per irrep */
    int *openpi;              /* number of singly-occupied MOs per irrep */
    int *virtpi;              /* number of vacant MOs per irrep */
    int *frozen_docc;         /* number of frozen doubly-occupied MOs per irrep */
    int *frozen_uocc;         /* number of frozen vacant MOs per irrep */
    int num_mo;               /* number of MOs */
    int ndocc;                /* number of doubly-occupied MOs */
    int nfrdocc;              /* number of "frozen" doubly occupied MOs */
    int nactdocc;             /* number of correlated doubly occupied MOs */
    int nsocc;                /* number of singly-occupied MOs */
    int nuocc;                /* number of vacant MOs */
    int nfruocc;              /* number of "frozen" vacant MOs */
    int nactuocc;             /* number of "active" vacant MOs */
    int num_moshells;         /* number of shells of MOs */
    int num_openmoshells;     /* number of shells of singly-occupied MOs */
} MOInfo_t;

