struct MOInfo {
    int nirreps;           /* no. of irreducible representations */
    int nmo;               /* no. of molecular orbitals */
    int noeints;           /* no. unique one-electron integrals (ex. fruocc) */
    int iopen;             /* 0=closed shell; >0=open shell */
    int *orbspi;           /* no. of MOs per irrep */
    int *clsdpi;           /* no. of closed-shells per irrep ex. frdocc */
    int *openpi;           /* no. of open-shells per irrep */
    int *uoccpi;           /* no. of unoccupied orbitals per irrep ex. fruocc */
    int *frdocc;           /* no. of frozen core orbitals per irrep */
    int *fruocc;           /* no. of frozen unoccupied orbitals per irrep */
    char **labels;         /* irrep labels */
    int nfzc;              /* total no. of frozen core orbitals */
    int nfzv;              /* total no. of frozen virtual orbitals */
    int nactive;           /* total no. of active orbitals */
    int *occ;              /* boolean array for active occ. orbs. */
    int *vir;              /* boolean array for active virt. orbs. */
    int *socc;             /* boolean array for active socc. orbs. */
    int *all_occ;          /* boolean array for occ. orbs. (in the full space) */
    int *all_vir;          /* boolean array for virt. orbs. (in the full space) */
    int *all_socc;         /* boolean array for socc. orbs. (in the full space) */
    int *frozen;           /* boolean array for frz. orbs (in the full space) */
    int *cc_occ;           /* QT->CC active occupied reordering array */
    int *cc_vir;           /* QT->CC active virtiual reordering array */
    int *cc_allocc;        /* QT->CC all occupied reordering array */
    int *cc_allvir;        /* QT->CC all virtual reordering array */
    int *qt_occ;           /* CC->QT active occupied reordering array */
    int *qt_vir;           /* CC->QT active virtiual reordering array */
    int *qt_allocc;        /* CC->QT all occupied reordering array */
    int *qt_allvir;        /* CC->QT all virtual reordering array */
    int *occ_sym;          /* CC active occupied index symmetry */
    int *vir_sym;          /* CC active virtual index symmetry */
    int *allocc_sym;       /* CC all occupied index symmetry */
    int *allvir_sym;       /* CC all virtual index symmetry */
    int *occpi;            /* no. of occupied orbs. (incl. open) per irrep */
    int *virtpi;           /* no. of virtual orbs. (incl. open) per irrep */
    int *all_occpi;        /* no. of occ. orbs. (incl. open and fzc) per irrep */
    int *all_virtpi;       /* no. of virt. orbs. (incl. open and fzc) per irrep */
    int *occ_off;          /* active occ. orbital offsets within each irrep */
    int *vir_off;         /* active virt. orbital offsets within each irrep */
    int *all_occ_off;      /* all occ. orbital offsets within each irrep */
    int *all_vir_off;     /* all virt. orbital offsets within each irrep */
    double enuc;           /* Nuclear repulsion energy */
    double efzc;           /* Frozen core energy */
    double eref;           /* The reference energy (computed here) */
};
