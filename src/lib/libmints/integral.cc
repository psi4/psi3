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

template <class T>
static void swap(T& x, T& y) {
    T tmp=x; x = y; y = tmp;
}

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

ShellCombinationsIterator IntegralFactory::shells_iterator()
{
    return ShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

IntegralsIterator ShellCombinationsIterator::integrals_iterator() 
{
    return IntegralsIterator(bs1_->shell(p()), bs2_->shell(q()), bs3_->shell(r()), bs4_->shell(s()));
}


void ShellCombinationsIterator::generate_combinations(const Ref<BasisSet> &bs1, const Ref<BasisSet> &bs2, const Ref<BasisSet> &bs3, const Ref<BasisSet> &bs4)
{
    // Assumes bs1 == bs2 == bs3 == bs4
    int usii, usjj, uskk, usll;
    int usi, usj, usk, usl;
    int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
    int num_unique_pk;
    
    for (usii=0; usii<bs1->nshell(); usii++) {
        for (usjj=0; usjj<=usii; usjj++) { 
            for (uskk=0; uskk<=usjj; uskk++) {
                for (usll=0; usll<=uskk; usll++) {
                    // Decide what shell quartets out of (ij|kl), (ik|jl), and (il|jk) are unique
                    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
                    if (usii == usjj && usii == uskk || usjj == uskk && usjj == usll)
                        num_unique_pk = 1;
                    else if (usii == uskk || usjj == usll) {
                        num_unique_pk = 2;
                        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
                    }
                    else if (usjj == uskk) {
                        num_unique_pk = 2;
                        usi_arr[1] = usii; usj_arr[1] = usll; usk_arr[1] = usjj; usl_arr[1] = uskk;
                    }
                    else if (usii == usjj || uskk == usll) {
                        num_unique_pk = 2;
                        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
                    }
                    else {
                        num_unique_pk = 3;
                        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
                        usi_arr[2] = usii; usj_arr[2] = usll; usk_arr[2] = usjj; usl_arr[2] = uskk;
                    }
                    
                    // For each num_unique_pk
                    for (int upk=0; upk < num_unique_pk; ++upk) {
                        usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];
                        
                        // Sort shells based on AM, saves ERI work doing permutation resorting.
                        if (bs1->shell(usi)->am(0) < bs2->shell(usj)->am(0)) {
                            swap(usi, usj);
                        }
                        if (bs3->shell(usk)->am(0) < bs4->shell(usl)->am(0)) {
                            swap(usk, usl);
                        }
                        if (bs1->shell(usi)->am(0) + bs2->shell(usj)->am(0) >
                            bs3->shell(usk)->am(0) + bs4->shell(usl)->am(0)) {
                            swap(usi, usk);
                            swap(usj, usl);
                        }
                        
                        ShellQuartet q;
                        q.P = usi; q.Q = usj; q.R = usk; q.S = usl;
                        unique_quartets_.push_back(q);
                    }
                }
            }
        }
    }    
}

void IntegralsIterator::generate_combinations(const Ref<GaussianShell> &usi, const Ref<GaussianShell> &usj, const Ref<GaussianShell> &usk, const Ref<GaussianShell> &usl)
{
    int ni =usi->nfunction(0);
    int nj =usj->nfunction(0);
    int nk =usk->nfunction(0);
    int nl =usl->nfunction(0);
    
    int fii = usi->function_index();
    int fij = usj->function_index();
    int fik = usk->function_index();
    int fil = usl->function_index();
    
    // For a given quartet save the individual integrals information.
    // Unfortunately we don't know if an individual integral exists until we compute it. :(
    if (usi == usj && usk == usl && usi == usk) { // (aa|aa) case
        int iimax = ni - 1;
        for (int ii=0; ii <= iimax; ++ii) {
            int jjmax = ii;
            for (int jj=0; jj <= jjmax; ++jj) {
                int kkmax = ii;
                for (int kk=0; kk <= kkmax; ++kk) {
                    int llmax = (kk==ii) ? jj : kk;
                    for (int ll=0; ll <= llmax; ++ll) {
                        // Save this integral
                        Integral integral;
                        integral.i = ii + fii;
                        integral.j = jj + fij;
                        integral.k = kk + fik;
                        integral.l = ll + fil;
                        integral.index = ll+nl*(kk+nk*(jj+nj*ii));
                        unique_integrals_.push_back(integral);
                    }
                }
            }
        }
    }
    else if (usi == usk && usj == usl) { // (ab|ab) case
        int iimax = ni-1;
        for (int ii=0; ii <= iimax; ++ii) {
            int jjmax = nj-1;
            for (int jj=0; jj <= jjmax; ++jj) {
                int kkmax = ii;
                for (int kk=0; kk <= kkmax; ++kk) {
                    int llmax = (kk==ii)? jj : nl - 1;
                    for (int ll=0; ll <= llmax; ++ll) {
                        Integral integral;
                        integral.i = ii + fii;
                        integral.j = jj + fij;
                        integral.k = kk + fik;
                        integral.l = ll + fil;
                        integral.index = ll+nl*(kk+nk*(jj+nj*ii));
                        
                        // Might need to swap indices
                        if (integral.i < integral.j) {
                            swap(integral.i, integral.j);
                            swap(integral.k, integral.l);
                        }
                        if (integral.i < integral.k) {
                            swap(integral.i, integral.k);
                            swap(integral.j, integral.l);
                        }
                        unique_integrals_.push_back(integral);
                    }
                }
            }
        }
    }
    else { // (ab|cd)
        int iimax = ni-1;
        int kkmax = nk-1;
        for (int ii=0; ii <= iimax; ++ii) {
            int jjmax = (usi == usj) ? ii : nj - 1;
            for (int jj=0; jj <= jjmax; ++jj) {
                for (int kk=0; kk <= kkmax; ++kk) {
                    int llmax = (usk == usl) ? kk : nl - 1;
                    for (int ll=0; ll <= llmax; ++ll) {
                        Integral integral;
                        integral.i = ii + fii;
                        integral.j = jj + fij;
                        integral.k = kk + fik;
                        integral.l = ll + fil;
                        integral.index = ll+nl*(kk+nk*(jj+nj*ii));
                        
                        // Might need to swap indices
                        if (integral.i < integral.j) {
                            swap(integral.i, integral.j);
                        }
                        if (integral.k < integral.l) {
                            swap(integral.k, integral.l);
                        }
                        if ((integral.i < integral.k) || (integral.i == integral.k && integral.j < integral.l)) {
                            swap(integral.i, integral.k);
                            swap(integral.j, integral.l);
                        }
                        unique_integrals_.push_back(integral);
                    }
                }
            }
        }
    }    
}
