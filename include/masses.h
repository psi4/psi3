/*
** The new PSI standard header file for atomic masses.
**
** Created by Matt Leininger, July 1995
**
** This file contains three arrays.  The first is a list of atomic masses
** for a number of isotopes.  The second is a list of isotope labels which
** correspond to the masses in the first array.  The mass labels given
** are the same as those in W. D. Allen's INTDER95.  The final array is a 
** list of atomic masses for the first several elements, with the most
** common isotopomer mass given.  This last array is most useful for
** converting an atomic number into an atomic mass. 
**
** Please be sure that any modifications to this file are reflected in
** _all three_ arrays if necessary.
**
** Modifications:
** 10/08/99 EFV - Changed an2masses[0] to 0.0000 
**                (ghost atom's weight should be 0)
**
** 08/16/02 EFV - Updated atomic_masses[] and an2masses for H-Ar using NIST database
**                of Atomic Weights and Isotopic Compositions (www.nist.gov)
**
** 08/23/02 EFV - Use lithium isotope 7 as the default (most abundant) 
**
*/

#ifndef _psi_include_masses_h_
#define _psi_include_masses_h_

#define LAST_MASS_INDEX (139)

#ifdef __cplusplus
extern "C" {
#endif

double atomic_masses[] = 
{ 1.0078250321, 1.0078250321, 2.0141017780, 3.0160492675, 2.0141017780, 3.0160492675, 
  4.0026032497, 3.0160293097, 4.0026032497,
  7.0160040, 6.0151223, 7.0160040, 
  9.0121821, 9.0121821,
 11.0093055, 10.0129370, 11.0093055,
 12.000000, 12.000000, 13.0033548378, 14.003241988,
 14.0030740052, 14.0030740052, 15.0001088984,
 15.9949146221, 15.9949146221, 16.99913150, 17.9991604,
 18.99840320, 18.99840320,
 19.9924401759, 19.9924401759, 20.99384674, 21.99138551,
 22.98976967, 22.98976967,
 23.98504190, 23.98504190, 24.98583702, 25.98259304,
 26.98153844, 26.98153844,
 27.9769265327, 27.9769265327, 28.97649472, 29.97377022,
 30.97376151, 30.97376151,
 31.97207069, 31.97207069, 32.97145850, 33.96786683, 35.96708088,
 34.96885271, 34.96885271, 36.96590260,
 39.962383123, 35.96754628, 37.9627322, 39.962383123,
 38.963708, 38.963708, 39.963999, 40.961825,
 39.962591, 39.962591, 41.958622, 42.958770, 43.955485, 45.953689, 47.952532,
 44.955914, 44.955914,
 47.947947, 45.952633, 46.951765, 47.947947, 48.947871, 49.944786,
 50.943963, 49.947161, 50.943963,
 51.940510, 49.946046, 51.940510, 52.940651, 53.938882,
 54.938046, 54.938046,
 55.934939, 53.939612, 55.934939, 56.935396, 57.933278,
 58.933198, 58.933198,
 57.935347, 57.935347, 59.930789, 60.931059, 61.928346, 63.927968,
 62.929599, 62.929599, 64.927792,
 63.929145, 63.929145, 65.926035, 66.927129, 67.924826, 69.925325,
 68.925581, 68.925581, 70.924701,
 73.921179, 69.924250, 71.922080, 72.923464, 73.921179, 75.921403,
 74.921596, 74.921596,
 79.916521, 73.922477, 75.919207, 76.919908, 77.917304, 79.916521, 81.916709,
 78.918336, 78.918336, 80.916290,
 83.911506, 77.920397, 79.916375, 81.913483, 82.914134, 83.911506, 85.910614} ;

char *mass_labels[] = 
{ "H", "H1", "H2", "H3", "D", "T",
  "HE", "HE3", "HE4",
  "LI", "LI6", "LI7",
  "BE", "BE9",
  "B", "B10", "B11",
  "C", "C12", "C13", "C14",
  "N", "N14", "N15",
  "O", "O16", "O17", "O18",
  "F", "F19",
  "NE", "NE20", "NE21", "NE22",
  "NA", "NA23",
  "MG", "MG24", "MG25", "MG26",
  "AL", "AL27",
  "SI", "SI28", "SI29", "SI30",
  "P", "P31",
  "S", "S32", "S33", "S34", "S36",
  "CL", "CL35", "CL37",
  "AR", "AR36", "AR38", "AR40",
  "K", "K39", "K40", "K41", 
  "CA", "CA40", "CA42", "CA43", "CA44", "CA46", "CA48", 
  "SC", "SC45",
  "TI", "TI46", "TI47", "TI48", "TI49", "TI50",
  "V", "V50", "V51",
  "CR", "CR50", "CR52", "CR53", "CR54",
  "MN", "MN55",
  "FE", "FE54", "FE56", "FE57", "FE58",
  "CO", "CO59",
  "NI", "NI58", "NI60", "NI61", "NI62", "NI64",
  "CU", "CU63", "CU65",
  "ZN", "ZN64", "ZN66", "ZN67", "ZN68", "ZN70",
  "GA", "GA69", "GA71",
  "GE", "GE70", "GE72", "GE73", "GE74", "GE76",
  "AS", "AS75",
  "SE", "SE74", "SE76", "SE77", "SE78", "SE80", "SE82",
  "BR", "BR79", "BR81",
  "KR", "KR78", "KR80", "KR82", "KR83", "KR84", "KR86"} ;

#define LAST_ATOMIC_INDEX (36)
char *atomic_labels[] = 
{ "X",
  "H", 
  "HE",
  "LI",
  "BE",
  "B", 
  "C", 
  "N", 
  "O", 
  "F", 
  "NE",
  "NA",
  "MG",
  "AL",
  "SI",
  "P", 
  "S", 
  "CL",
  "AR",
  "K", 
  "CA",
  "SC",
  "TI",
  "V", 
  "CR",
  "MN",
  "FE",
  "CO",
  "NI",
  "CU",
  "ZN",
  "GA",
  "GE",
  "AS",
  "SE",
  "BR",
  "KR" };


double an2masses[] =
{ 0.0000000, 1.0078250321, 4.0026032497,  6.0151223, 9.0121821, 11.0093055, 12.000000,
  14.0030740052, 15.9949146221, 18.99840320, 19.9924401759, 22.98976967,  23.98504190, 
  26.98153844, 27.9769265327, 30.97376151, 31.97207069, 34.96885271, 39.962383123,
  38.963708, 39.962591, 44.955914, 47.947947, 50.943963, 51.940510, 
  54.938046, 55.934939, 58.933198, 57.935347, 62.929599, 63.929145,
  68.925581, 73.921179, 74.921596, 79.916521, 78.918336, 83.911506
} ; 

#ifdef __cplusplus
}
#endif

#endif /* header guard */
