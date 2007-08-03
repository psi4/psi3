/*! \file 
    \ingroup (INTDER)
    \brief Enter brief description of file here 
*/
#ifndef __cartesian_h
#define __cartesian_h

class Cartesian
{
  double x, y, z;

public:
  Cartesian(double x1, double y1, double z1)
    { x = x1; y = y1; z = z1; }
  Cartesian()
    { x = 0.0; y = 0.0; z = 0.0; }

  double Distance(Cartesian& B);
  double Distance(double bx, double by, double bz);

  double& getX()
    { return x; }
  double& getY()
    { return y; }
  double& getZ()
    { return z; }
};

#endif

