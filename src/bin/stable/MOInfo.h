struct MOInfo {
    int nirreps;           /* no. of irreducible representations */
    int nmo;               /* no. of molecular orbitals */
    int iopen;             /* 0=closed shell; >0=open shell */
    int *orbspi;           /* no. of MOs per irrep */
    int *clsdpi;           /* no. of closed-shells per irrep  */
    int *openpi;           /* no. of open-shells per irrep */
    int *uoccpi;           /* no. of unoccupied orbitals per irrep  */
    char **labels;         /* irrep labels */
    int *occ;              /* boolean array for active occ. orbs. */
    int *vir;              /* boolean array for active virt. orbs. */
    int *socc;             /* boolean array for active socc. orbs. */
    int *cc_occ;           /* QT->CC active occupied reordering array */
    int *cc_vir;           /* QT->CC active virtiual reordering array */
    int *qt_occ;           /* CC->QT active occupied reordering array */
    int *qt_vir;           /* CC->QT active virtiual reordering array */
    int *occ_sym;          /* CC active occupied index symmetry */
    int *vir_sym;          /* CC active virtual index symmetry */
    int *occpi;            /* no. of occupied orbs. (incl. open) per irrep */
    int *virtpi;           /* no. of virtual orbs. (incl. open) per irrep */
    int *occ_off;          /* occ. orbital offsets within each irrep */
    int *vir_off;          /* virt. orbital offsets within each irrep */
    int *frdocc;
    int *fruocc;
    int nfzc;
    int nfzv;
};
