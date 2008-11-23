#ifndef _psi_src_lib_libmints_onebody_h_
#define _psi_src_lib_libmints_onebody_h_

/*!
    \file libmints/onebody.h
    \ingroup MINTS
*/

#include <libmints/ref.h>
#include <libmints/matrix.h>

namespace psi {
    
class IntegralFactory;
class BasisSet;
class GaussianShell;

/// Basis class for all one-electron integrals.
class OneBodyInt
{
protected:
    IntegralFactory *integral_;
    Ref<BasisSet> bs1_;
    Ref<BasisSet> bs2_;
    
    double *buffer_;
    unsigned int count_;
    int deriv_;
    int natom_;
    
    OneBodyInt(IntegralFactory *integral, const Ref<BasisSet>& bs1, const Ref<BasisSet>& bs2=0, int deriv=0);
    
public:
    virtual ~OneBodyInt();
    
    /// Basis set on center one.
    Ref<BasisSet> basis();
    /// Basis set on center one.
    Ref<BasisSet> basis1();
    /// Basis set on center two.
    Ref<BasisSet> basis2();
    
    /// Buffer where the integrals are placed.
    const double *buffer() const;
    
    /// Compute the integrals between basis function in the given shell pair.
    virtual void compute_shell(int, int) = 0;
    
    /// Computes all integrals and stores them in result
    void compute(RefMatrix& result);
    
    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(RefMatrixArray& result);
    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(RefSimpleMatrixArray& result);
    
    /// Does the method provide first derivatives?
    virtual bool has_deriv1() { return false; }
    
    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(RefMatrixArray& result);
    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(RefSimpleMatrixArray& result);
    
    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv1(int, int);
    
    /// Integral object that created me.
    IntegralFactory *integral() const { return integral_; }
    
    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();
    
    /// Returns a clone of this object. By default throws an exception.
    virtual Ref<OneBodyInt> clone();
    
    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(Ref<GaussianShell> & , Ref<GaussianShell> &, int nchunk=1);
    
    /// Transform Cartesian integrals to spherical harmonic ones.
    /// Reads from buffer_ and stores results back in buffer_.
    virtual void spherical_transform(Ref<GaussianShell> & , Ref<GaussianShell> &);
    
    void do_transform(Ref<GaussianShell> &, Ref<GaussianShell> &, int);
    
    /// Accumulates results into a RefMatrix
    void so_transform(RefMatrix       &result, int, int, int ichunk=0);
    /// Accumulates results into a RefSimpleMatrix
    void so_transform(RefSimpleMatrix &result, int, int, int ichunk=0);
};

}

#endif
