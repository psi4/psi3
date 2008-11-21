#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <psifiles.h>

#include "basisset.h"
#include "integral.h"
#include "symmetry.h"

using namespace psi;

BasisSet::BasisSet(Ref<Chkpt> &chkpt)
{
    max_nprimitives_ = 0;
    nshells_      = chkpt->rd_nshell();
    nprimitives_  = chkpt->rd_nprim();
    nao_          = chkpt->rd_nao();
    max_stability_index_ = 0;
    
    // Psi3 only allows either all Cartesian or all Spherical harmonic
    nbf_          = chkpt->rd_nso();
    puream_       = chkpt->rd_puream() ? true : false;
    max_am_       = chkpt->rd_max_am();
    uso2ao_       = chkpt->rd_usotao();
    
    // Allocate memory for the shells
    shells_       = new Ref<GaussianShell>[nshells_];
    
    // Initialize the shells
    initialize_shells(chkpt);
}

BasisSet::~BasisSet()
{
    Chkpt::free(shell_first_basis_function_);
    Chkpt::free(shell_first_ao_);
    Chkpt::free(shell_center_);
    Chkpt::free(uso2ao_);
}

void BasisSet::initialize_shells(Ref<Chkpt> &chkpt)
{
    // Retrieve angular momentum of each shell (1=s, 2=p, ...)
    int *shell_am = chkpt->rd_stype();
    
    // Retrieve number of primitives per shell
    int *shell_num_prims = chkpt->rd_snumg();
    
    // Retrieve exponents of primitive Gaussians
    double *exponents = chkpt->rd_exps();
    
    // Retrieve coefficients of primitive Gaussian
    double **ccoeffs = chkpt->rd_contr_full();
    
    // Retrieve pointer to first primitive in shell
    int *shell_fprim = chkpt->rd_sprim();
    
    // Retrieve pointer to first basis function in shell
    shell_first_basis_function_ = chkpt->rd_sloc_new();
    
    // Retrieve pointer to first AO in shell
    shell_first_ao_ = chkpt->rd_sloc();
    
    // Retrieve location of shells (which atom it's centered on)
    shell_center_ = chkpt->rd_snuc();
    
    // Initialize molecule, retrieves number of center and geometry
    molecule_ = new Molecule;
    molecule_->init_with_chkpt(chkpt);
    
    // Initialize SphericalTransform
    for (int i=0; i<=max_am_; ++i) {
        sphericaltransforms_.push_back(SphericalTransform(i));
    }
    
    // Initialize SOTransform
    sotransform_ = new SOTransform;
    sotransform_->init(nshells_);
    
    int *so2symblk = new int[nbf_];
    int *so2index  = new int[nbf_];
    int *sopi = chkpt->rd_sopi();
    int nirreps = chkpt->rd_nirreps();
    
    // Create so2symblk and so2index
    int ij = 0; int offset = 0;
    for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<sopi[h]; ++i) {
            so2symblk[ij] = h;
            so2index[ij] = ij-offset;
            
            ij++;
        }
        offset += sopi[h];
    }
    
    // Currently all basis sets are treated as segmented contractions
    // even though GaussianShell is generalized (well not really).
    int ncontr = 1;
    int ao_start = 0;
    int puream_start = 0;
    int *sym_transform = new int[nirreps];
    
    // We need access to symmetry information found in the checkpoint file.
    Symmetry symmetry(chkpt);
    
    for (int i=0; i<nshells_; ++i) {
        int *am = new int[ncontr];
        am[0] = shell_am[i] - 1;
        int fprim = shell_fprim[i] - 1;
        int nprims = shell_num_prims[i];
        Vector3 center = molecule_->xyz(shell_center_[i] - 1);
        double **cc = new double*[nprims];
        for (int p=0; p<nprims; ++p) {
            cc[p] = new double[ncontr];
            cc[p][0] = ccoeffs[fprim+p][am[0]];
        }
        
        // Construct a new shell. GaussianShell copies the data to new memory
        shells_[i] = new GaussianShell(ncontr, nprims, &(exponents[fprim]), am, 
            puream_ ? GaussianShell::Pure : GaussianShell::Cartesian, cc, shell_center_[i]-1, center,
            puream_start);
            
        if (nprims > max_nprimitives_)
            max_nprimitives_ = nprims;
            
        for (int p=0; p<nprims; p++) {
            delete[] cc[p];
        }
        delete[] cc;
        
        // OK, for a given number of AO functions in a shell INT_NCART(am)
        // beginning at column ao_start go through all rows finding where this
        // AO function contributes to an SO.
        for (int ao = 0; ao < INT_NCART(am[0]); ++ao) {
            int aooffset = ao_start + ao;
            for (int so = 0; so < nbf_; ++so) {
                if (fabs(uso2ao_[so][aooffset]) >= 1.0e-14)
                    sotransform_->add_transform(i, so2symblk[so], so2index[so], uso2ao_[so][aooffset], ao, so);
            }
        }
        
        // Set the symmetry transform vector of the shell
        for (int j=0; j<nirreps; ++j) {
            sym_transform[j] = symmetry.trans_vec(i, j);
        }
        // The shell will copy the elements into a new vector.
        shells_[i]->set_sym_transform(nirreps, sym_transform);
    
        // Compute index of the stabilizer for the shell
        int count = 1;
        for (int g=1; g<nirreps; ++g) {
            if (i == symmetry.trans_vec(i, g)-1)
                count++;
        }
        int stab_index = nirreps / count;
        if (max_stability_index_ < stab_index)
            max_stability_index_ = stab_index;
        
        // Shift the ao starting index over to the next shell
        ao_start += INT_NCART(am[0]);
        puream_start += INT_NFUNC(puream_, am[0]);
    }
    
    delete[] sym_transform;
    delete[] so2symblk;
    delete[] so2index;
    free(sopi);
    free_block(ccoeffs);
    free(exponents);
    free(shell_am);
    free(shell_num_prims);
    free(shell_fprim);
}

void BasisSet::print(FILE *out) const
{
    fprintf(out, "  Basis Set\n");
    fprintf(out, "    Number of shells: %d\n", nshell());
    fprintf(out, "    Number of basis function: %d\n", nbf());
    fprintf(out, "    Number of Cartesian functions: %d\n", nao());
    fprintf(out, "    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
    fprintf(out, "    Max angular momentum: %d\n\n", max_am());
    
    fprintf(out, "    Shells:\n\n");
    for (int s=0; s<nshell(); ++s)
        shells_[s]->print(out);
}

Ref<GaussianShell>& BasisSet::shell(int si) const
{
    #ifdef DEBUG
    assert(si < nshell());
    #endif
    return shells_[si];
}
