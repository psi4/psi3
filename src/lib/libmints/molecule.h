#ifndef _psi_src_lib_libmints_molecule_h_
#define _psi_src_lib_libmints_molecule_h_

#include <vector>
#include <string>
#include <cstdio>

#include <libmints/ref.h>
#include <libmints/vector3.h>
#include <libmints/vector.h>

#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

extern FILE *outfile;

namespace psi {
    
class Molecule
{
public:
    typedef struct atom_info {
        double x, y, z;
        int Z;
        double charge;
        double mass;
        std::string label;
    };
    
protected:
    /// Number of atoms.
    int natoms_;
    /// Atom info vector
    std::vector<atom_info> atoms_;
    /// Symmetry information about the molecule
    int nirreps_;
    /// Zero it out
    void clear();
    
public:
    Molecule();
    virtual ~Molecule();
    
    /// Pull information from a chkpt object created from psio
    void init_with_chkpt(Ref<psi::PSIO> &psio);
    /// Pull information from the chkpt object passed
    void init_with_chkpt(Ref<psi::Chkpt> &chkpt);
    
    /// Add an atom to the molecule
    void add_atom(int Z, double x, double y, double z,
                  const char * = 0, double mass = 0.0,
                  int have_charge = 0, double charge = 0.0);

    /// Number of atoms
    int natom() const { return natoms_; }
    /// Nuclear charge of atom
    int Z(int atom) const { return atoms_[atom].Z; }
    /// Return reference to atom_info struct for atom
    const atom_info &r(int atom) const { return atoms_[atom]; }
    /// Return copy of atom_info for atom
    atom_info r(int atom) { return atoms_[atom]; }
    /// Returns a Vector3 with x, y, z position of atom
    const Vector3 xyz(int atom) const { return Vector3(atoms_[atom].x, atoms_[atom].y, atoms_[atom].z); }
    /// Returns mass atom atom
    double mass(int atom) const;
    /// Returns label of atom
    const std::string label(int atom) const;
    /// Returns charge of atom
    double charge(int atom) const { return atoms_[atom].charge; }
    
    /// Tests to see of an atom is at the passed position with a given tolerance
    int atom_at_position(double *, double tol = 0.05) const;
    
    /// Computes center of mass of molecule (does not translate molecule)
    Vector3 center_of_mass() const;
    /// Computes nuclear repulsion energy
    double nuclear_repulsion_energy();
    /// Returns the nuclear contribution to the dipole moment
    SimpleVector nuclear_dipole_contribution();
    /// Returns the nuclear contribution to the quadrupole moment
    SimpleVector nuclear_quadrupole_contribution();
    
    /// Translates molecule by r
    void translate(const Vector3& r);
    /// Moves molecule to center of mass
    void move_to_com();
    
    /// Returns the number of irreps
    int nirrep() const { return nirreps_; }
    /// Sets the number of irreps
    void nirrep(int nirreps) { nirreps_ = nirreps; }
    
    /// Print the molecule
    void print(FILE *out = outfile);
};

}

#endif