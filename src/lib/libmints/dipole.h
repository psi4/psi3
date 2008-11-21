#ifndef _psi_src_lib_libmints_dipole_h_
#define _psi_src_lib_libmints_dipole_h_

#include <libmints/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

namespace psi {
/// Computes dipole integrals.
/// Use an IntegralFactory to create this object.
class DipoleInt : public OneBodyInt
{
    ObaraSaikaTwoCenterRecursion overlap_recur_;
    
    void compute_pair(Ref<GaussianShell> &, Ref<GaussianShell> &);
    void compute_pair_deriv1(Ref<GaussianShell> &, Ref<GaussianShell> &);
        
public:
    DipoleInt(IntegralFactory*, Ref<BasisSet> &, Ref<BasisSet> &, int deriv=0);
    virtual ~DipoleInt();
    
    void compute_shell(int, int);
    void compute_shell_deriv1(int, int);
    
    void compute(RefSimpleMatrixArray& result);
    void compute_deriv1(RefSimpleMatrixArray& result);
    
    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

}

#endif
