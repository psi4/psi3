struct moinfo {
  int nmo;               /* no. of MOs */
  int nso;               /* no. of symmetry adapted MOs */
  int nirreps;           /* no. of irreducible representations */
  char **irreplabels;    /* irrep labels */
  int *mopi;             /* all MOs per irrep */
  int *doccpi;           /* all doubly occupied MOs per irrep */
  int *soccpi;           /* all singly occupied MOs per irrep */
  int nfzdocc;           /* no. of frozen occupied MOs */
  int nfzvirt;           /* no. of frozen virtual MOs */
  int nactmo;            /* no. of active MOs */
  int *fzdoccpi;         /* frozen occupied MOs per irrep */
  int *fzvirtpi;         /* frozen virtual MOs per irrep */
  
  int *occpi;            /* occupied MOs per irrep */
  int *virpi;            /* virtual MOs per irrep */
  int *occ_sym;          /* occupied MOs symmetry */
  int *vir_sym;          /* virtual MOs symmetry */
  int *occ_off;          /* occupied orbital offsets within each irrep */ 
  int *vir_off;          /* virtual orbital offsets within each irrep */
  int *qt_occ;           /* CC->QT active occupied reordering array */
  int *qt_vir;           /* CC->QT active virtual reordering array */
  
  int *aoccpi;           /* alpha occupied orbitals per irrep */
  int *boccpi;           /* beta occupied orbitals per irrep */
  int *avirpi;           /* alpha virtual orbitals per irrep */
  int *bvirpi;           /* beta virtual orbitals per irrep */
  int *aocc_sym;         /* alpha occupied MOs symmetry */
  int *bocc_sym;         /* beta occupied MOs symmetry */
  int *avir_sym;         /* alpha virtual MOs symmetry */
  int *bvir_sym;         /* beta virtual MOs symmetry */
  int *aocc_off;         /* alpha occupied orbital offsets within each irrep */
  int *bocc_off;         /* beta occupied orbital offsets within each irrep */
  int *avir_off;         /* alpha virtual orbital offsets within each irrep */
  int *bvir_off;         /* beta virtual orbital offsets within each irrep */
  int *qt_aocc;          /* CC->QT alpha occupied reordering array */
  int *qt_bocc;          /* CC->QT beta occupied reordering array */
  int *qt_avir;          /* CC->QT alpha virtual reordering array*/
  int *qt_bvir;          /* CC->QT beta virtual reordering array */
  
  int *ioff;             /* ioff array */
  double Enuc;           /* Nuclear repulsion energy */
  double Escf;           /* SCF energy */
  double Emp2;           /* MP2 energy */
};
