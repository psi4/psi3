/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/

/*! \defgroup CCDENSITY ccdensity: Computes the Coupled-Cluster Density */

namespace psi { namespace ccdensity {

struct Frozen {
    int nfzc;           /* total no. of frozen core orbitals */
    int nfzv;           /* total no. of frozen virtual orbitals */
    int *occ_sym;       /* occupied index symmetry */
    int *aocc_sym;      /* alpha occupied index symmetry */
    int *bocc_sym;      /* beta occupied index symmetry */
    int *vir_sym;       /* virtual index symmetry */
    int *avir_sym;      /* alpha virtual index symmetry */
    int *bvir_sym;      /* beta virtual index symmetry */
    int *occpi;         /* no. of occ. orbs. (incl. open) per irrep */
    int *aoccpi;        /* no. of alpha occ. orbs. (incl. open) per irrep */
    int *boccpi;        /* no. of beta occ. orbs. (incl. open) per irrep */
    int *virtpi;        /* no. of virt. orbs. (incl. open) per irrep */
    int *avirtpi;       /* no. of alpha virt. orbs. (incl. open) per irrep */
    int *bvirtpi;       /* no. of beta virt. orbs. (incl. open) per irrep */
    int *occ_off;       /* occupied orbital offsets within each irrep */
    int *aocc_off;      /* occupied alpha orbital offsets within each irrep */
    int *bocc_off;      /* occupied beta orbital offsets within each irrep */
    int *vir_off;       /* virtual orbital offsets within each irrep */
    int *avir_off;      /* virtual alpha orbital offsets within each irrep */
    int *bvir_off;      /* virtual beta orbital offsets within each irrep */
    int *allcc_occ;     /* QT->CC occupied reordering array */
    int *allcc_vir;     /* QT->CC virtiual reordering array */
    int *qt_occ;        /* CC->QT occupied reordering array */
    int *qt_vir;        /* CC->QT virtiual reordering array */
    int *cc_occ;        /* QT->CC active occupied reordering array */
    int *cc_vir;        /* QT->CC active virtiual reordering array */
    int *occ;           /* boolean array for occ. orbs. */
    int *vir;           /* boolean array for virt. orbs. */
    int *socc;          /* boolean array for socc. orbs. */
};

}} // namespace psi::ccdensity
