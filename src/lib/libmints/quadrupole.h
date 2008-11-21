#ifndef _psi_src_lib_libmints_quadrupole_h_
#define _psi_src_lib_libmints_quadrupole_h_

#include <libmints/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

namespace psi {
    
/// Computes quadrupole integrals. At last check this may not be working.
/// Use an IntegralFactory to create this object.
class QuadrupoleInt : public OneBodyInt
{
    ObaraSaikaTwoCenterRecursion overlap_recur_;
    
    void compute_pair(Ref<GaussianShell> &, Ref<GaussianShell> &);
    
public:
    QuadrupoleInt(IntegralFactory*, Ref<BasisSet> &, Ref<BasisSet> &);
    virtual ~QuadrupoleInt();
    
    void compute_shell(int, int);
    
    /// Computes all quadrupole integrals (Qxx, Qxy, Qxz, Qyy, Qyz, Qzz) result must be an array of enough
    /// size to contain it.
    void compute(RefMatrixArray& result);
    void compute(RefSimpleMatrixArray& result);
    
    virtual void spherical_transform(Ref<GaussianShell> & , Ref<GaussianShell> &);
};

}

#endif
