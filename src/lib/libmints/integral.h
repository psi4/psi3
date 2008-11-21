#ifndef _psi_src_lib_libmints_integral_h_
#define _psi_src_lib_libmints_integral_h_

#include <libmints/ref.h>
#include <vector>

#define INT_NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)
#define INT_NPURE(am) (2*(am)+1)
#define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))
#define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))

namespace psi {
    
class BasisSet;
class OneBodyInt;
class TwoBodyInt;
class Symmetry;

class SphericalTransformComponent
{
protected:
    int a_, b_, c_;
    int cartindex_, pureindex_;
    
    double coef_;
    
public:
    /// Returns the exponent of x.
    int a() const { return a_; }
    /// Returns the exponent of y.
    int b() const { return b_; }
    /// Returns the exponent of z.
    int c() const { return c_; }
    /// Returns the index of the Cartesian basis function
    int cartindex() const { return cartindex_; }
    /// Returns the index of the spherical harmonic basis function
    int pureindex() const { return pureindex_; }
    /// Returns the coefficient of this component of the transformation
    double coef() const { return coef_; }
    
    void init(int a, int b, int c, double coef, int cartindex, int pureindex);
};

class SphericalTransform
{
protected:
    std::vector<SphericalTransformComponent> components_;
    int l_; // The angular momentum this transform is for.
    
    SphericalTransform();
public:
    SphericalTransform(int l);
    virtual ~SphericalTransform() {};
        
    /// Returns the Cartesian basis function index of component i
    int cartindex(int i) const { return components_[i].cartindex(); }
    /// Returns the spherical harmonic basis index of component i
    int pureindex(int i) const { return components_[i].pureindex(); }
    /// Returns the transformation coefficient of component i
    double coef(int i) const { return components_[i].coef(); }
    /// Returns the Cartesian basis function's x exponent of component i
    int a(int i) const { return components_[i].a(); }
    /// Returns the Cartesian basis function's y exponent of component i
    int b(int i) const { return components_[i].b(); }
    /// Returns the Cartesian basis function's z exponent of component i
    int c(int i) const { return components_[i].c(); }
    /// Returns the number of components in the transformation
    int n() const { return components_.size(); }
    /// Returns the angular momentum
    int l() const { return l_; }
};

class SphericalTransformIter
{
private:
    SphericalTransform *trans_;
    int i_;
    
public:
    SphericalTransformIter(SphericalTransform* trans) { trans_ = trans; i_ = 0; }
    
    void first() { i_ = 0; }
    void next()  { i_++;   }
    bool is_done() { return i_ < trans_->n() ? true : false; }
    
    /// Returns the Cartesian basis function index of component i
    int cartindex() const { return trans_->cartindex(i_); }
    /// Returns the spherical harmonic basis index of component i
    int pureindex() const { return trans_->pureindex(i_); }
    /// Returns the transformation coefficient of component i
    double coef()   const { return trans_->coef(i_); }
    /// Returns the Cartesian basis function's x exponent of component i
    int a()         const { return trans_->a(i_); }
    /// Returns the Cartesian basis function's y exponent of component i
    int b()         const { return trans_->b(i_); }
    /// Returns the Cartesian basis function's z exponent of component i
    int c()         const { return trans_->c(i_); }
};

class IntegralFactory
{
protected:
    /// Center 1 basis set
    Ref<BasisSet> bs1_;
    /// Center 2 basis set
    Ref<BasisSet> bs2_;
    /// Center 3 basis set
    Ref<BasisSet> bs3_;
    /// Center 4 basis set
    Ref<BasisSet> bs4_;
    
    /// Provides ability to transform to and from sphericals (d=0, f=1, g=2)
    std::vector<SphericalTransform> spherical_transforms_;
    
public:
    /** Initialize IntegralFactory object given a GaussianBasisSet for each center. */
    IntegralFactory(const Ref<BasisSet> &bs1, const Ref<BasisSet> &bs2,
                    const Ref<BasisSet> &bs3, const Ref<BasisSet> &bs4);
    
    virtual ~IntegralFactory();
    
    /// Set the basis set for each center.
    virtual void set_basis(const Ref<BasisSet> &bs1, const Ref<BasisSet> &bs2 = 0,
        const Ref<BasisSet> &bs3 = 0, const Ref<BasisSet> &bs4 = 0);
        
    /// Returns an OneBodyInt that computes the overlap integral.
    virtual Ref<OneBodyInt> overlap(int deriv=0);
    
    /// Returns an OneBodyInt that computes the kinetic energy integral.
    virtual Ref<OneBodyInt> kinetic(int deriv=0);
    
    /// Returns an OneBodyInt that computes the nuclear attraction integral.
    virtual Ref<OneBodyInt> potential(int deriv=0);

    /// Returns an OneBodyInt that computes the dipole integral.
    virtual Ref<OneBodyInt> dipole(int deriv=0);
    
    /// Returns an OneBodyInt that computes the quadrupole integral.
    virtual Ref<OneBodyInt> quadrupole();
    
    /// Returns an ERI integral object
    virtual Ref<TwoBodyInt> eri();

    /// Returns an ERI integral derivative object
    virtual Ref<TwoBodyInt> eri_deriv();
    
    /// Initializes spherical harmonic transformations
    virtual void init_spherical_harmonics(int max_am);
    
    // Return spherical transform object for am
    SphericalTransform* spherical_transform(int am) { return &(spherical_transforms_[am]); }
};

}

#endif