struct MOInfo {
    int nirreps;        /* no. of irreducible representations */
    int nmo;            /* no. of molecular orbitals */
    int iopen;          /* 0=closed shell; >0=open shell */
    int *orbspi;        /* no. of MOs per irrep */
    int *clsdpi;        /* no. of closed-shells per irrep excl. frdocc */
    int *openpi;        /* no. of open-shells per irrep */
    int *uoccpi;        /* no. of unoccupied orbitals per irrep excl. fruocc */
    int *frdocc;        /* no. of frozen core orbitals per irrep */
    int *fruocc;        /* no. of frozen unoccupied orbitals per irrep */
    char **labels;      /* irrep labels */
    int nfzc;           /* total no. of frozen core orbitals */
    int nfzv;           /* total no. of frozen virtual orbitals */
    int nclsd;          /* total no. of closd shells excl. frdocc */
    int nopen;          /* total no. of open shells  */
    int nuocc;          /* total no. of unoccupied shells excl. fruocc */
    int *occ_sym;       /* active occupied index symmetry */
    int *vir_sym;       /* active virtual index symmetry */
    int *occpi;         /* no. of active occ. orbs. (incl. open) per irrep */
    int *virtpi;        /* no. of active virt. orbs. (incl. open) per irrep */
    int *occ_off;       /* occupied orbital offsets within each irrep */
    int *vir_off;       /* virtual orbital offsets within each irrep */
    int *qt_occ;        /* CC->QT active occupied reordering array */
    int *qt_vir;        /* CC->QT active virtiual reordering array */
    double enuc;        /* Nuclear repulsion energy */
    double escf;        /* SCF energy from file30 */
    double eref;        /* Reference energy */
    double ecc;         /* CC energy from ccenergy */
    double **opdm;      /* Onepdm in the full (fzc+clsd+socc+uocc) space */
    double **I;         /* Lagrangian matrix in the full space */
};
