#ifndef __displacements_h
#define __displacements_h

#include "molecule.h"
#include <vector>

class Displacements
{
  std::vector<Molecule> vectorMolecules;

public:
  bool loadFromOptKing();
  bool loadFromCheckPoint();
  bool loadFromInput();

  void printGeometries();
  void moveToCenterOfMass();

  void atom_num(char*, double*);
  void addDisplacement(Molecule& mol);
  void useMasses(double *mass);
  Molecule* displacement(int disp);
};


#endif

