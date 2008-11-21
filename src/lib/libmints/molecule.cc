#include <cmath>
#include <cstdio>

#include <libmints/molecule.h>

#include <masses.h>
#include <physconst.h>

using namespace std;
using namespace psi;

Molecule::Molecule():
    natoms_(0), nirreps_(0)
{
    
}

Molecule::~Molecule()
{
    clear();
}

void Molecule::clear()
{
    natoms_ = 0;
    nirreps_ = 0;
    atoms_.empty();
}

void Molecule::add_atom(int Z, double x, double y, double z,
                        const char *label, double mass,
                        int have_charge, double charge)
{
    atom_info info;
    
    info.x = x;
    info.y = y;
    info.z = z;
    info.Z = Z;
    info.charge = charge;
    info.label = label;
    info.mass = mass;
    
    natoms_++;
    atoms_.push_back(info);
}

double Molecule::mass(int atom) const
{
    if (atoms_[atom].mass != 0.0)
        return atoms_[atom].mass;
        
    return an2masses[atoms_[atom].Z];
}

const std::string Molecule::label(int atom) const
{
    return atoms_[atom].label; 
}

int Molecule::atom_at_position(double *xyz, double tol) const
{
    Vector3 b(xyz);
    for (int i=0; i < natom(); ++i) {
        Vector3 a(atoms_[i].x, atoms_[i].y, atoms_[i].z);
        if (b.distance(a) < tol)
            return i;
    }
    return -1;
}

Vector3 Molecule::center_of_mass() const
{
    Vector3 ret;
    double total_m;
    
    ret = 0.0;
    total_m = 0.0;
    
    for (int i=0; i<natom(); ++i) {
        double m = mass(i);
        ret += m * xyz(i);
        total_m += m;
    }
    
    ret *= 1.0/total_m;
    
    return ret;
}

double Molecule::nuclear_repulsion_energy()
{
    double e=0.0;
    
    for (int i=1; i<natom(); ++i) {
        for (int j=0; j<i; ++j) {
            e += charge(i) * charge(j) / (xyz(i).distance(xyz(j)));
        }
    }
    
    return e;
}

void Molecule::translate(const Vector3& r)
{
    for (int i=0; i<natom(); ++i) {
        atoms_[i].x += r[0];
        atoms_[i].y += r[1];
        atoms_[i].z += r[2];
    }
}

void Molecule::move_to_com()
{
    Vector3 com = -center_of_mass();
    translate(com);
}

void Molecule::init_with_chkpt(Ref<PSIO> &psio)
{
    // User sent a psio object. Create a chkpt object based on it.
    Ref<Chkpt> chkpt(new Chkpt(psio.pointer(), PSIO_OPEN_OLD));
    init_with_chkpt(chkpt);
}

void Molecule::init_with_chkpt(Ref<Chkpt> &chkpt)
{
    int atoms = chkpt->rd_natom();
    double *zvals = chkpt->rd_zvals();
    double **geom = chkpt->rd_geom();
    
    for (int i=0; i<atoms; ++i) {
        add_atom((int)zvals[i], geom[i][0], geom[i][1], geom[i][2], atomic_labels[(int)zvals[i]], an2masses[(int)zvals[i]]);
    }
    
    nirreps_ = chkpt->rd_nirreps();
    
    Chkpt::free(zvals);
    Chkpt::free(geom);
}

void Molecule::print(FILE *out)
{
    if (natom()) {
        fprintf(out,"       Center              X                  Y                   Z\n");
        fprintf(out,"    ------------   -----------------  -----------------  -----------------\n");
        
        for(int i = 0; i < natom(); ++i){
            Vector3 geom = xyz(i);
            fprintf(out, "    %12s ",label(i).c_str()); fflush(out);
            for(int j = 0; j < 3; j++)
                fprintf(out, "  %17.12f", geom[j]);
            fprintf(out,"\n");
        }
        fprintf(out,"\n");
        fflush(out);
    }
}

SimpleVector Molecule::nuclear_dipole_contribution()
{
    SimpleVector ret(3);
    
    for(int i=0; i<natom(); ++i) {
        Vector3 geom = xyz(i);
        ret[0] += Z(i) * geom[0];
        ret[1] += Z(i) * geom[1];
        ret[2] += Z(i) * geom[2];
    }
        
    return ret;
}

SimpleVector Molecule::nuclear_quadrupole_contribution()
{
    SimpleVector ret(6);
    double xx, xy, xz, yy, yz, zz;
    
    xx = xy = xz = yy = yz = zz = 0.0;
    
    for (int i=0; i<natom(); ++i) {
        Vector3 geom = xyz(i);
        ret[0] += Z(i) * geom[0] * geom[0]; // xx
        ret[1] += Z(i) * geom[0] * geom[1]; // xy
        ret[2] += Z(i) * geom[0] * geom[2]; // xz
        ret[3] += Z(i) * geom[1] * geom[1]; // yy
        ret[4] += Z(i) * geom[1] * geom[2]; // yz
        ret[5] += Z(i) * geom[2] * geom[2]; // zz
    }
            
    return ret;
}
