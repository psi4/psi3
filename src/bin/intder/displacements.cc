/*! \file 
    \ingroup (INTDER)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include "displacements.h"
#include "params.h"

extern "C" {
  #include <libchkpt/chkpt.h>
  #include <libciomr/libciomr.h>
  #include <libipv1/ip_lib.h>
  #include <libpsio/psio.h>
  #include <psifiles.h>
  #include <masses.h>
  #include <physconst.h>

  extern FILE* outfile;
  extern Params gParams;
  extern char *psi_file_prefix;
};

void Displacements::addDisplacement(Molecule& mol)
{
  vectorMolecules.push_back(mol);
}

Molecule* Displacements::displacement(int disp)
{
  return &(vectorMolecules[disp]);
}

void Displacements::printGeometries()
{
  int index;
  
  for (index = 0; index < vectorMolecules.size(); index++)
    {
      if (index == 0)
	fprintf(outfile, "Disp #%d (Reference Geometry)\n", index);
      else
	fprintf(outfile, "Disp #%d\n", index);
      vectorMolecules[index].printGeometry();
    }
}

void Displacements::moveToCenterOfMass()
{
  int index;
  
  for (index = 0; index < vectorMolecules.size(); index++)
    vectorMolecules[index].moveToCenterOfMass();
}

void Displacements::useMasses(double *mass)
{
  int index;

  for (index = 0; index < vectorMolecules.size(); index++)
    vectorMolecules[index].useMasses(mass);
}

