struct Frozen {
    int nfzc;           /* total no. of frozen core orbitals */
    int nfzv;           /* total no. of frozen virtual orbitals */
    int *occ_sym;       /* occupied index symmetry */
    int *vir_sym;       /* virtual index symmetry */
    int *occpi;         /* no. of occ. orbs. (incl. open) per irrep */
    int *virtpi;        /* no. of virt. orbs. (incl. open) per irrep */
    int *occ_off;       /* occupied orbital offsets within each irrep */
    int *vir_off;       /* virtual orbital offsets within each irrep */
    int *allcc_occ;     /* QT->CC occupied reordering array */
    int *allcc_vir;     /* QT->CC virtiual reordering array */
    int *qt_occ;        /* CC->QT occupied reordering array */
    int *qt_vir;        /* CC->QT virtiual reordering array */
    int *cc_occ;        /* QT->CC active occupied reordering array */
    int *cc_vir;        /* QT->CC active virtiual reordering array */
    int *occ;              /* boolean array for occ. orbs. */
    int *vir;              /* boolean array for virt. orbs. */
    int *socc;             /* boolean array for socc. orbs. */
};
