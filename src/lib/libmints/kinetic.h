#ifndef _psi_src_lib_libmints_kinetic_h_
#define _psi_src_lib_libmints_kinetic_h_

#include <libmints/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/integral.h>

namespace psi {
    
class KineticInt : public OneBodyInt
{
    ObaraSaikaTwoCenterRecursion overlap_recur_;
    
    void compute_pair(Ref<GaussianShell> &, Ref<GaussianShell> &);
    void compute_pair_deriv1(Ref<GaussianShell> &, Ref<GaussianShell> &);
    
public:
    KineticInt(IntegralFactory*, Ref<BasisSet> &, Ref<BasisSet> &, int deriv=0);
    virtual ~KineticInt();
    
    void compute_shell(int, int);
    void compute_shell_deriv1(int, int);
    
    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
};

}

#endif
