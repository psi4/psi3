/*! \file atom.h
    \ingroup (INTDER)
    \brief Enter brief description of file here 
*/
#ifndef __atom_h
#define __atom_h

#include "cartesian.h"

class Atom : public Cartesian
{
public:
  int atomicNumber;
  double atomicWeight;
  Atom(int an, double aw, double ax, double ay, double az);
  Atom(int an, double aw, Cartesian& A);
  Atom(int an, double aw);
  Atom();

  int getAtomicNumber()
    { return atomicNumber; }
  void setAtomicNumber(int an)
    { atomicNumber = an; }
  double getAtomicWeight()
    { return atomicWeight; }
  void setAtomicWeight(double aw)
    { atomicWeight = aw; }
};

#endif

