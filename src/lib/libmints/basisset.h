#ifndef _psi_src_lib_libmints_basisset_h_
#define _psi_src_lib_libmints_basisset_h_

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libmints/ref.h>
#include <libmints/molecule.h>
#include <libmints/gshell.h>
#include <libmints/sobasis.h>
#include <libmints/integral.h>

extern FILE *outfile;

namespace psi {
class BasisSet
{
    int nprimitives_;
    int nshells_;
    int nao_;
    int nbf_;
    int max_am_;
    int max_nprimitives_;
    int *shell_first_basis_function_;
    int *shell_first_ao_;
    int *shell_center_;
    int max_stability_index_;
    double **uso2ao_;
    
    bool puream_;
    
    Ref<Ref<GaussianShell>, SimpleReferenceCount, StandardArrayPolicy> shells_;
    Ref<Molecule> molecule_;
    Ref<SOTransform> sotransform_;
    std::vector<SphericalTransform> sphericaltransforms_;
    
    // No default constructor
    BasisSet();
    // No assignment
    BasisSet& operator=(const BasisSet&);
    
    void initialize_shells(Ref<psi::Chkpt> &chkpt);
    
public:
    /// Constructor, reads in the basis set from the checkpoint file
    BasisSet(Ref<psi::Chkpt> &chkpt);
    /// Copy constructor, currently errors if used
    BasisSet(const BasisSet&);
    /// Destructor
    ~BasisSet();
    
    /// Total number of primitives
    int nprimitive() const             { return nprimitives_; }
    /// Maximum number of primitives in a shell
    int max_nprimitive() const         { return max_nprimitives_; }
    /// Number of shells
    int nshell() const                 { return nshells_;     }
    /// Number of atomic orbitals
    int nao() const                    { return nao_;         }
    /// Number of basis functions
    int nbf() const                    { return nbf_;         }
    /// Maximum angular momentum
    int max_am() const                 { return max_am_;      }
    /// Spherical harmonics?
    bool has_puream() const            { return puream_;      }
    /// Molecule this basis is for
    Ref<Molecule> molecule() const     { return molecule_;    }
    /// Maximum stabilizer index
    int max_stability_index() const    { return max_stability_index_; }
    /// Given a shell what is its first AO function
    int shell_to_function(int i) const { return shell_first_ao_[i]; }
    int shell_to_basis_function(int i) const { return shell_first_basis_function_[i]; }
    
    /// Return the si'th Gaussian shell
    Ref<GaussianShell>& shell(int si) const;
    
    /// Returns i'th shell's transform
    SOTransformShell* so_transform(int i) { return sotransform_->aoshell(i); }
    
    /// Returns the transformation object for a given angular momentum. Used in ERIs.
    SphericalTransform* spherical_transform(int am) { return &sphericaltransforms_[am]; }
    
    /// Print the basis set
    void print(FILE *out = outfile) const;
};

}

#endif
