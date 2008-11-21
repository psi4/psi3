#ifndef _psi_src_lib_libmints_twobody_h
#define _psi_src_lib_libmints_twobody_h

#include <libmints/ref.h>
#include <libmints/matrix.h>

namespace psi {
    
class IntegralFactory;
class BasisSet;
class GaussianShell;

class TwoBodyInt
{
protected:
    IntegralFactory *integral_;
    Ref<BasisSet> bs1_;
    Ref<BasisSet> bs2_;
    Ref<BasisSet> bs3_;
    Ref<BasisSet> bs4_;
    
    /// Buffer to hold the final integrals.
    double *target_;
    /// Buffer to hold the transformation intermediates.
    double *tformbuf_;
    /// Buffer to hold the initially computed integrals.
    double *source_;
    /// Maximum number of unique quartets needed to compute a set of SO's
    int max_unique_quartets_;
    
    void permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24);
    void permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    
    TwoBodyInt(IntegralFactory *integral,
               const Ref<BasisSet>& bs1,
               const Ref<BasisSet>& bs2,
               const Ref<BasisSet>& bs3,
               const Ref<BasisSet>& bs4);
               
public:
    virtual ~TwoBodyInt();
    
    /// Basis set on center one
    Ref<BasisSet> basis();
    /// Basis set on center one
    Ref<BasisSet> basis1();
    /// Basis set on center two
    Ref<BasisSet> basis2();
    /// Basis set on center three
    Ref<BasisSet> basis3();
    /// Basis set on center four
    Ref<BasisSet> basis4();

    /// Buffer where the integrals are placed
    const double *buffer() const { return target_; };
    
    /// Compute the integrals
    virtual void compute_shell(int, int, int, int) = 0;
    
    /// Integral object that created me.
    IntegralFactory *integral() const { return integral_; }
    
    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(Ref<GaussianShell> &, Ref<GaussianShell> &, Ref<GaussianShell> &, Ref<GaussianShell> &, int nchunk=1);
        
    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();
    
    /// Returns a clone of this object. By default throws an exception
    virtual Ref<TwoBodyInt> clone();
    
    /// Results go back to buffer_
    void pure_transform(int, int, int, int, int ichunk=0);
};

}

#endif
