#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/overlap.h>
#include <libmints/kinetic.h>
#include <libmints/potential.h>
#include <libmints/integral.h>
#include <libmints/dipole.h>
#include <libmints/quadrupole.h>
#include <libmints/symmetry.h>
#include <libmints/eri.h>

using namespace psi;

/** Initialize IntegralFactory object given a GaussianBasisSet for each center. */
IntegralFactory::IntegralFactory(const Ref<BasisSet> &bs1, const Ref<BasisSet> &bs2,
                const Ref<BasisSet> &bs3, const Ref<BasisSet> &bs4)
{
    set_basis(bs1, bs2, bs3, bs4);
}

IntegralFactory::~IntegralFactory()
{
    
}

void IntegralFactory::set_basis(const Ref<BasisSet> &bs1, const Ref<BasisSet> &bs2,
    const Ref<BasisSet> &bs3, const Ref<BasisSet> &bs4)
{
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;
    
    // Find the max am
    Ref<BasisSet> max12 = bs1_->max_am() > bs2_->max_am() ? bs1_ : bs2_;
    Ref<BasisSet> max34 = bs3_->max_am() > bs4_->max_am() ? bs3_ : bs4_;
    Ref<BasisSet> max1234 = max12->max_am() > max34->max_am() ? max12 : max34;
    
    init_spherical_harmonics(max1234->max_am());
}

Ref<OneBodyInt> IntegralFactory::overlap(int deriv)
{
    return new OverlapInt((IntegralFactory*)this, bs1_, bs2_, deriv);
}

Ref<OneBodyInt> IntegralFactory::kinetic(int deriv)
{
    return new KineticInt((IntegralFactory*)this, bs1_, bs2_, deriv);
}

Ref<OneBodyInt> IntegralFactory::potential(int deriv)
{
    return new PotentialInt((IntegralFactory*)this, bs1_, bs2_, deriv);
}

Ref<OneBodyInt> IntegralFactory::dipole(int deriv)
{
    return new DipoleInt((IntegralFactory*)this, bs1_, bs2_, deriv);
}

Ref<OneBodyInt> IntegralFactory::quadrupole()
{
    return new QuadrupoleInt((IntegralFactory*)this, bs1_, bs2_);
}

Ref<TwoBodyInt> IntegralFactory::eri(int deriv)
{
    return new ERI((IntegralFactory*)this, bs1_, bs2_, bs3_, bs4_, deriv);
}

void IntegralFactory::init_spherical_harmonics(int max_am)
{
    for (int i=0; i<=max_am; ++i)
        spherical_transforms_.push_back(SphericalTransform(i));
}
