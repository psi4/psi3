#include <stdlib.h>
#include <stdio.h>
#include "cartesian.h"
#include "atom.h"
#include "molecule.h"
#include "displacements.h"
#include "params.h"
#include "intco.h"
#include "globals.h"
//#include "bmat.h"
//#include "transform.h"

extern "C" {
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>

char *gprgid(void);
};

// Global variables
Displacements gDisplacements;
Params gParams;
InternalCoordinates gIntCo;
//BMat gBMat;
//Transform gTransform;

void intro(void);

int main(int argc, char** argv)
{
  double *e12 = NULL;
  double *e23 = NULL;
  double *e13 = NULL;
  double distance = 0.0;
  double angle = 0.0;
  int index;
  int i,j;

  psi_start(argc-1, argv+1, 0);
  ip_cwk_add(gprgid());
  psio_init();
  tstart(outfile);
  chkpt_init(PSIO_OPEN_OLD);

  intro();

  gParams.parseInputFile();
  gParams.printParams();

  //Reading geometry
  gDisplacements.loadFromCheckPoint();
  //gDisplacements.loadFromOptKing();
  //gDisplacements.loadFromInput();
  
  if((gParams.analysisType != 0) && (gParams.analysisType != 5) || gParams.propertyDimension > 0 || gParams.generateDisp < 0 
     || gParams.n_rcom || gParams.transformType < 0) {
    gParams.checkMasses();
    gParams.checkIsotopes(); // isotopes will overwrite masses, but isotopes returns if !ISOTOPES
  }
  
  gParams.writeMasses();

  fprintf(outfile, "\nNuclear Cartesian coordinates in angstroms:\n");
  gDisplacements.printGeometries();

  //Need to save center of mass transformation variables? Need to be able to save non-center of mass
  //transformed coordinates

  if((gParams.transformType > 3 || gParams.numtest == 3) && gParams.ndisp == 0) {
    gDisplacements.moveToCenterOfMass();
    fprintf(outfile, "\nCoordinates shifted to center of mass\n"); 
  }
  
  fflush(outfile);
 
  /*
    fprintf(outfile, "\nB Matrix for Internal Coordinates:\n");
    gBMat.disp = 0;
    gBMat.init();
    gBMat.make_BMat();
    
    //Once symmcoords are implemented, this logic will need to be changed
    if(!gParams.symmCoords) {
    for(i = 0; i < gParams.nintco; i++) {
    gBMat.SVectSave[i] = gBMat.SVectArray[i];
    for(j = 0; j < gParams.ncartesians; j++)
    gBMat.BMatSave[i][j] = gBMat.BMatrix[i][j];
    }
    }
    
    fprintf(outfile,"\n\tReprinting saved simple internal coordinates with bmat s-values (angstroms or degrees)\n");
    for(i = 0; i < gParams.nintco; i++) 
    gIntCo.printSingleIntCo(i);
    
    if(gParams.print_array[0] > 2)
    gBMat.print_BMat();
    
    if(gParams.analysisType >= 0) {  //else goes to line 600, aka frequency analysis
    gBMat.invert_BMat();
    
    if(gParams.derlvl == 0) {
    fprintf(outfile, "\nderlvl = 0, no transformation requested");
    exit(0);
    }
    if(gParams.derlvl == 1 && gParams.atEquilibrium == 0) {
    fprintf(outfile, "\nAt equilibrium, all gradients = 0)");
    exit(0);
    } 
    
    //Intder now makes blocks of memory to determine dimensions of 2/3/4-D matrices. We shouldn't have to do this 
    
    //This doesn't work yet
    //if(gParams.transformType == 2)
    //readIntCoDerivatives();
    
    if(gParams.derlvl + gParams.atEquilibrium >= 3) {
    SecondDerivative secondDer;
    if(gParams.numtest > 1)
    secondDer.NumTest();
    
    secondDer.i_to_c();
    
    if(gParams.numtest == 1)
    secondDer.SRTest1();
    }
    
    if(gParams.derlvl + gParams.atEquilibrium >= 4) {
    ThirdDerivative thirdDer;
    
    if(abs(gParams.numtest == 2))
    thirdDer.NumTest();
    
    thirdDer.i_to_c();
    if(gParams.numtest == 2)
    thirdDer.SRTest2();
    
    }
    
    if(gParams.derlvl + gParams.atEquilibrium >= 5) {
    FourthDerivative fourthDer;
    
    fourthDer.NumTest();
    fourthDer.i_to_c();
    }
    }  
    
    if(gParams.stop != 1) {
    
    if(gParams.propertyDimension > 0) 
    //  XKinput();
    //if(gParams.transformType < 0) 
    //linear_transform();
    
    }
    */
  
  chkpt_close();
  tstop(outfile);
  psio_done();
  psi_stop();
  return 0;
}

char *gprgid(void)
{
  char *prgid = ":INTDER";
  return (prgid);
}

void intro(void)
{
    fprintf(outfile,
            "\n\t------------------------------------------------------\n");
    fprintf(outfile,
              "\t                    INTDER2C5                         \n");
    fprintf(outfile,
              "\t   Developed by Wesley D. Allen and Co-workers.       \n");
    fprintf(outfile,
              "\t   PSI 3.2 C++ Implementation by Nathan DeYonker,     \n");
    fprintf(outfile,
              "\t   Justin Turney, and Rollin King.                    \n");
    fprintf(outfile,
              "\t   Performs various vibrational analyses and higher-  \n");
    fprintf(outfile,
              "\t   order nonlinear transformations among force field  \n");
    fprintf(outfile,
              "\t   representations.                                   \n");
    fprintf(outfile,
              "\t------------------------------------------------------\n");
}