#include "atom.h"
#include "molecule.h"
#include "params.h"

extern "C" {
extern double atomic_masses[];
extern char *mass_labels[];
extern char *atomic_labels[];
extern double an2masses[];
extern FILE* outfile;

extern Params gParams;
};

Molecule::Molecule()
{
}

Molecule::~Molecule()
{
}

void Molecule::moveToCenterOfMass()
{
  double xm = 0.0;
  double ym = 0.0;
  double zm = 0.0;
  double mt = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  int index = 0;

  for (index = 0; index < vectorAtoms.size(); index++)  {
    xm += vectorAtoms[index].getAtomicWeight() * vectorAtoms[index].getX();
    ym += vectorAtoms[index].getAtomicWeight() * vectorAtoms[index].getY();
    zm += vectorAtoms[index].getAtomicWeight() * vectorAtoms[index].getZ();
    mt += vectorAtoms[index].getAtomicWeight(); 
  }

  cx = xm / mt;
  cy = ym / mt;
  cz = zm / mt;

  for (index = 0; index < vectorAtoms.size(); index++)
  {
    vectorAtoms[index].getX() -= cx;
    vectorAtoms[index].getY() -= cy;
    vectorAtoms[index].getZ() -= cz;
  }
}

void Molecule::addAtom(int an, double aq, double ax, double ay, double az)
{
  vectorAtoms.push_back(Atom(an, aq, ax, ay, az));
}

Atom* Molecule::atom(int an)
{
  return &(vectorAtoms[an]);
}

void Molecule::printGeometry()
{
  int index = 0;

  for (index = 0; index < vectorAtoms.size(); index++)
  {
    fprintf(outfile, "Atom #%d %s(%15.10lf): %15.10lf %15.10lf %15.10lf\n", index+1,
              atomic_labels[vectorAtoms[index].getAtomicNumber()],
              vectorAtoms[index].getAtomicWeight(),
              vectorAtoms[index].getX(),
              vectorAtoms[index].getY(),
              vectorAtoms[index].getZ());
  }
}

void Molecule::useMasses(double *mass)
{
  int index = 0;

  for (index = 0; index < gParams.natom; index++) {
      vectorAtoms[index].setAtomicWeight(mass[index]);
  }
}

