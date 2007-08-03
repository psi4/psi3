/*! \file 
    \ingroup (INTDER)
    \brief Enter brief description of file here 
*/
#ifndef __molecule_h
#define __molecule_h

#include "atom.h"
#include <vector>

class Molecule
{
  std::vector<Atom> vectorAtoms;

public:
  Molecule();
  ~Molecule();

  void addAtom(int an, double aw, double ax, double ay, double az);
  void useMasses(double *mass);
  void moveToCenterOfMass();

  void printGeometry();
  Atom* atom(int an);
};

#endif

