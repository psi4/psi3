struct MOInfo {
    int nirreps;        /* no. of irreducible representations */
    int nmo;            /* no. of molecular orbitals */
    int iopen;          /* 0=closed shell; >0=open shell */
    int phase;          /* Boolean for consistency of orbital phases */
    int *orbspi;        /* no. of MOs per irrep */
    int *clsdpi;        /* no. of closed-shells per irrep excl. frdocc */
    int *openpi;        /* no. of open-shells per irrep */
    int *uoccpi;        /* no. of unoccupied orbitals per irrep excl. fruocc */
    int *frdocc;        /* no. of frozen core orbitals per irrep */
    int *fruocc;        /* no. of frozen unoccupied orbitals per irrep */
    char **labels;      /* irrep labels */
    int *occ_sym;       /* relative occupied index symmetry */
    int *vir_sym;       /* relative virtual index symmetry */
    int *occpi;         /* no. of occupied orbs. (incl. open) per irrep */
    int *virtpi;        /* no. of virtual orbs. (incl. open) per irrep */
    int *occ_off;       /* occupied orbital offsets within each irrep */
    int *vir_off;       /* virtual orbital offsets within each irrep */
};
