#include "blas.h"
#include "mrcc.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::perturbative_cbs()
{

  fprintf(outfile,"\n\n  Computing perturbative corrections for frozen-virtual MOs");

  blas->solve("t2_1[oo][vv]{u}   = <[oo]:[vv]> / d2[oo][vv]{u}");
  blas->solve("t2_1[oO][vV]{u}   = <[oo]|[vv]> / d2[oO][vV]{u}");
  blas->solve("t2_1[OO][VV]{u}   = <[oo]:[vv]> / d2[OO][VV]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     t2_1[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VV]{u} . <[oo]:[vv]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vv = blas->get_scalar("ECCSD",0);

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_1[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  blas->solve("t2_1[oo][vf]{u}   = <[oo]:[vf]> / d2[oo][vf]{u}");
  blas->solve("t2_1[oO][vF]{u}   = <[oo]|[vf]> / d2[oO][vF]{u}");
  blas->solve("t2_1[OO][VF]{u}   = <[oo]:[vf]> / d2[OO][VF]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_1[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf = blas->get_scalar("ECCSD",0);

  blas->solve("t2_1[oo][fv]{u}   = <[oo]:[fv]> / d2[oo][fv]{u}");
  blas->solve("t2_1[oO][fV]{u}   = <[oo]|[fv]> / d2[oO][fV]{u}");
  blas->solve("t2_1[OO][FV]{u}   = <[oo]:[fv]> / d2[OO][FV]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][fv]{u} . <[oo]:[fv]>");
  blas->solve("Eabab{u} =     t2_1[oO][fV]{u} . <[oo]|[fv]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][FV]{u} . <[oo]:[fv]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_fv = blas->get_scalar("ECCSD",0);


  blas->solve("t2_1[oo][ff]{u}   = <[oo]:[ff]> / d2[oo][ff]{u}");
  blas->solve("t2_1[oO][fF]{u}   = <[oo]|[ff]> / d2[oO][fF]{u}");
  blas->solve("t2_1[OO][FF]{u}   = <[oo]:[ff]> / d2[OO][FF]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_1[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff = blas->get_scalar("ECCSD",0);

  fprintf(outfile,"\n\n        CBS(virt-fvir)           = %20.12f",E_vf);
  fprintf(outfile,"\n        CBS(fvir-virt)           = %20.12f",E_fv);
  fprintf(outfile,"\n        CBS(fvir-fvir)           = %20.12f",E_ff);
  fprintf(outfile,"\n\n      * CBS total correction     = %20.12f",E_ff + E_vf + E_fv);


  // Sort amplitudes
  blas->solve("t2_1[ov][ov]{u} = #1324# t2_1[oo][vv]{u}");
  blas->solve("t2_1[ov][OV]{u} = #1324# t2_1[oO][vV]{u}");
  blas->solve("t2_1[oV][Ov]{u} = #1342# t2_1[oO][vV]{u}");
  blas->solve("t2_1[OV][ov]{u} = #3142# t2_1[oO][vV]{u}");
  blas->solve("t2_1[OV][OV]{u} = #1324# t2_1[OO][VV]{u}");

  blas->solve("t2_1[ov][of]{u} = #1324# t2_1[oo][vf]{u}");
  blas->solve("t2_1[ov][OF]{u} = #1324# t2_1[oO][vF]{u}");
  blas->solve("t2_1[OV][OF]{u} = #1324# t2_1[OO][VF]{u}");
  blas->solve("t2_1[OV][of]{u} = #3142# t2_1[oO][fV]{u}");

  blas->solve("t2_1[of][ov]{u} = #1324# t2_1[oo][fv]{u}");
  blas->solve("t2_1[of][OV]{u} = #1324# t2_1[oO][fV]{u}");
  blas->solve("t2_1[OF][ov]{u} = #3142# t2_1[oO][vF]{u}");
  blas->solve("t2_1[OF][OV]{u} = #1324# t2_1[OO][FV]{u}");

  blas->solve("t2_1[of][of]{u} = #1324# t2_1[oo][ff]{u}");
  blas->solve("t2_1[of][OF]{u} = #1324# t2_1[oO][fF]{u}");
  blas->solve("t2_1[OF][of]{u} = #3142# t2_1[oO][fF]{u}");
  blas->solve("t2_1[OF][OF]{u} = #1324# t2_1[OO][FF]{u}");

  blas->solve("t2_1[oF][Ov]{u} = #1342# t2_1[oO][vF]{u}");
  blas->solve("t2_1[oF][Of]{u} = #1342# t2_1[oO][fF]{u}");
  blas->solve("t2_1[oV][Of]{u} = #1342# t2_1[oO][fV]{u}");

  // OOVV Terms
  // HHHH Ladder
  blas->solve("t2_eqns[oo][vv]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[oo][vv]{u}");

  blas->solve("t2_eqns[oO][vV]{u}  =     <[oo]|[oo]> 1@1 t2_1[oO][vV]{u}");

  blas->solve("t2_eqns[OO][VV]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[OO][VV]{u}");

  // PPPP Ladder
  blas->solve("t2_eqns[oo][vv]{u} += 1/2 t2_1[oo][vv]{u} 2@2 <[vv]:[vv]>");
  blas->solve("t2_eqns[oo][vv]{u} += 1/2 t2_1[oo][vf]{u} 2@2 <[vv]:[vf]>");
  blas->solve("t2_eqns[oo][vv]{u} += 1/2 t2_1[oo][fv]{u} 2@2 <[vv]:[fv]>");
  blas->solve("t2_eqns[oo][vv]{u} += 1/2 t2_1[oo][ff]{u} 2@2 <[vv]:[ff]>");

  blas->solve("t2_eqns[oO][vV]{u} +=  t2_1[oO][vV]{u} 2@2 <[vv]|[vv]>");
  blas->solve("t2_eqns[oO][vV]{u} +=  t2_1[oO][vF]{u} 2@2 <[vv]|[vf]>");
  blas->solve("t2_eqns[oO][vV]{u} +=  t2_1[oO][fV]{u} 2@2 <[vv]|[fv]>");
  blas->solve("t2_eqns[oO][vV]{u} +=  t2_1[oO][fF]{u} 2@2 <[vv]|[ff]>");

  blas->solve("t2_eqns[OO][VV]{u} += 1/2 t2_1[OO][VV]{u} 2@2 <[vv]:[vv]>");
  blas->solve("t2_eqns[OO][VV]{u} += 1/2 t2_1[OO][VF]{u} 2@2 <[vv]:[vf]>");
  blas->solve("t2_eqns[OO][VV]{u} += 1/2 t2_1[OO][FV]{u} 2@2 <[vv]:[fv]>");
  blas->solve("t2_eqns[OO][VV]{u} += 1/2 t2_1[OO][FF]{u} 2@2 <[vv]:[ff]>");

  // HPHP Ring
  blas->solve("t2_eqns[oo][vv]{u} += #1342#   t2_1[ov][ov]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #1432# - t2_1[ov][ov]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #2341# - t2_1[ov][ov]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #2431#   t2_1[ov][ov]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #1342#   t2_1[ov][OV]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #1432# - t2_1[ov][OV]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #2341# - t2_1[ov][OV]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oo][vv]{u} += #2431#   t2_1[ov][OV]{u} 2@2 ([vo]|[ov])");

  blas->solve("t2_eqns[oo][vv]{u} += #1342#   t2_1[ov][of]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #1432# - t2_1[ov][of]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #2341# - t2_1[ov][of]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #2431#   t2_1[ov][of]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #1342#   t2_1[ov][OF]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #1432# - t2_1[ov][OF]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #2341# - t2_1[ov][OF]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oo][vv]{u} += #2431#   t2_1[ov][OF]{u} 2@2 ([vo]|[of])");

  blas->solve("t2_eqns[oO][vV]{u} += #1342#   t2_1[ov][ov]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oO][vV]{u} += #1342#   t2_1[ov][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oO][vV]{u} += #1423# - t2_1[oV][Ov]{u} 2@2 <[ov]|[ov]>");
  blas->solve("t2_eqns[oO][vV]{u} += #2314# - t2_1[oV][Ov]{u} 1@2 <[ov]|[ov]>");
  blas->solve("t2_eqns[oO][vV]{u} += #2431#   t2_1[OV][OV]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oO][vV]{u} += #2431#   t2_1[OV][ov]{u} 2@2 ([vo]:[ov])");

  blas->solve("t2_eqns[oO][vV]{u} += #1342#   t2_1[ov][of]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oO][vV]{u} += #1342#   t2_1[ov][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oO][vV]{u} += #1423# - t2_1[oV][Of]{u} 2@2 <[ov]|[of]>");
  blas->solve("t2_eqns[oO][vV]{u} += #2314# - t2_1[oF][Ov]{u} 1@2 <[ov]|[of]>");
  blas->solve("t2_eqns[oO][vV]{u} += #2431#   t2_1[OV][OF]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oO][vV]{u} += #2431#   t2_1[OV][of]{u} 2@2 ([vo]:[of])");

  blas->solve("t2_eqns[OO][VV]{u} += #1342#   t2_1[OV][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #1432# - t2_1[OV][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #2341# - t2_1[OV][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #2431#   t2_1[OV][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #1342#   t2_1[OV][ov]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #1432# - t2_1[OV][ov]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #2341# - t2_1[OV][ov]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[OO][VV]{u} += #2431#   t2_1[OV][ov]{u} 2@2 ([vo]|[ov])");

  blas->solve("t2_eqns[OO][VV]{u} += #1342#   t2_1[OV][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #1432# - t2_1[OV][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #2341# - t2_1[OV][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #2431#   t2_1[OV][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #1342#   t2_1[OV][of]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #1432# - t2_1[OV][of]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #2341# - t2_1[OV][of]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[OO][VV]{u} += #2431#   t2_1[OV][of]{u} 2@2 ([vo]|[of])");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][vv]{u}   = t2_eqns[oo][vv]{u} / d2[oo][vv]{u}");
  blas->solve("t2_2[oO][vV]{u}   = t2_eqns[oO][vV]{u} / d2[oO][vV]{u}");
  blas->solve("t2_2[OO][VV]{u}   = t2_eqns[OO][VV]{u} / d2[OO][VV]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     t2_2[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][VV]{u} . <[oo]:[vv]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vv_3 = blas->get_scalar("ECCSD",0);

  fprintf(outfile,"\n\n      * CBS total correction (2)  = %20.12f",E_vv);
  fprintf(outfile,"\n\n      * CBS total correction (3)  = %20.12f",E_vv_3);
  fprintf(outfile,"\n\n      * CBS total correction (2+3)= %20.12f",E_vv + E_vv_3);

  // OOVF Terms
  // HHHH Ladder
  blas->solve("t2_eqns[oo][vf]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[oo][vf]{u}");

  blas->solve("t2_eqns[oO][vF]{u}  =     <[oo]|[oo]> 1@1 t2_1[oO][vF]{u}");

  blas->solve("t2_eqns[OO][VF]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[OO][VF]{u}");

  // PPPP Ladder
  blas->solve("t2_eqns[oo][vf]{u} += 1/2 t2_1[oo][vv]{u} 2@2 <[vf]:[vv]>");
  blas->solve("t2_eqns[oo][vf]{u} += 1/2 t2_1[oo][vf]{u} 2@2 <[vf]:[vf]>");
  blas->solve("t2_eqns[oo][vf]{u} += 1/2 t2_1[oo][fv]{u} 2@2 <[vf]:[fv]>");
  blas->solve("t2_eqns[oo][vf]{u} += 1/2 t2_1[oo][ff]{u} 2@2 <[vf]:[ff]>");

  blas->solve("t2_eqns[oO][vF]{u} +=  t2_1[oO][vV]{u} 2@2 <[vf]|[vv]>");
  blas->solve("t2_eqns[oO][vF]{u} +=  t2_1[oO][vF]{u} 2@2 <[vf]|[vf]>");
  blas->solve("t2_eqns[oO][vF]{u} +=  t2_1[oO][fV]{u} 2@2 <[vf]|[fv]>");
  blas->solve("t2_eqns[oO][vF]{u} +=  t2_1[oO][fF]{u} 2@2 <[vf]|[ff]>");

  blas->solve("t2_eqns[OO][VF]{u} += 1/2 t2_1[OO][VV]{u} 2@2 <[vf]:[vv]>");
  blas->solve("t2_eqns[OO][VF]{u} += 1/2 t2_1[OO][VF]{u} 2@2 <[vf]:[vf]>");
  blas->solve("t2_eqns[OO][VF]{u} += 1/2 t2_1[OO][FV]{u} 2@2 <[vf]:[fv]>");
  blas->solve("t2_eqns[OO][VF]{u} += 1/2 t2_1[OO][FF]{u} 2@2 <[vf]:[ff]>");

  // HPHP Ring
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2_1[ov][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #1432# - t2_1[of][ov]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2_1[ov][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2431#   t2_1[of][ov]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2_1[ov][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #1432# - t2_1[of][OV]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2_1[ov][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2431#   t2_1[of][OV]{u} 2@2 ([vo]|[ov])");

  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2_1[of][ov]{u} 1@2 ([fo]:[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #1432# - t2_1[of][of]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2_1[of][ov]{u} 1@2 ([fo]:[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #2431#   t2_1[of][of]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2_1[ov][OF]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #1432# - t2_1[of][OF]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2_1[ov][OF]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oo][vf]{u} += #2431#   t2_1[of][OF]{u} 2@2 ([vo]|[of])");

  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2_1[ov][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2_1[ov][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #1423# - t2_1[oF][Ov]{u} 2@2 <[ov]|[ov]>");
  blas->solve("t2_eqns[oO][vF]{u} += #2314# - t2_1[oV][Ov]{u} 1@2 <[of]|[ov]>");
  blas->solve("t2_eqns[oO][vF]{u} += #2431#   t2_1[OF][OV]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #2431#   t2_1[OF][ov]{u} 2@2 ([vo]:[ov])");

  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2_1[of][ov]{u} 1@2 ([fo]|[of])");
  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2_1[ov][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[oO][vF]{u} += #1423# - t2_1[oF][Of]{u} 2@1 <[of]|[ov]>");
  blas->solve("t2_eqns[oO][vF]{u} += #2314# - t2_1[oF][Ov]{u} 1@2 <[of]|[of]>");
  blas->solve("t2_eqns[oO][vF]{u} += #2431#   t2_1[OF][OF]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[oO][vF]{u} += #2431#   t2_1[OF][of]{u} 2@2 ([vo]:[of])");

  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2_1[OV][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #1432# - t2_1[OF][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2_1[OV][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2431#   t2_1[OF][OV]{u} 2@2 ([vo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2_1[OV][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #1432# - t2_1[OF][ov]{u} 2@2 ([vo]|[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2_1[OV][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2431#   t2_1[OF][ov]{u} 2@2 ([vo]|[ov])");

  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2_1[OV][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #1432# - t2_1[OF][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2_1[OV][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #2431#   t2_1[OF][OF]{u} 2@2 ([vo]:[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2_1[OV][of]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #1432# - t2_1[OF][of]{u} 2@2 ([vo]|[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2_1[OV][of]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[OO][VF]{u} += #2431#   t2_1[OF][of]{u} 2@2 ([vo]|[of])");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][vf]{u}   = t2_eqns[oo][vf]{u} / d2[oo][vf]{u}");
  blas->solve("t2_2[oO][vF]{u}   = t2_eqns[oO][vF]{u} / d2[oO][vF]{u}");
  blas->solve("t2_2[OO][VF]{u}   = t2_eqns[OO][VF]{u} / d2[OO][VF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_2[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf_3 = blas->get_scalar("ECCSD",0);


  // OOFF Terms
  // HHHH Ladder
  blas->solve("t2_eqns[oo][ff]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[oo][ff]{u}");

  blas->solve("t2_eqns[oO][fF]{u}  =     <[oo]|[oo]> 1@1 t2_1[oO][fF]{u}");

  blas->solve("t2_eqns[OO][FF]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[OO][FF]{u}");

  // PPPP Ladder
  blas->solve("t2_eqns[oo][ff]{u} += 1/2 t2_1[oo][vv]{u} 2@2 <[ff]:[vv]>");
  blas->solve("t2_eqns[oo][ff]{u} += 1/2 t2_1[oo][vf]{u} 2@2 <[ff]:[vf]>");
  blas->solve("t2_eqns[oo][ff]{u} += 1/2 t2_1[oo][fv]{u} 2@2 <[ff]:[fv]>");
  blas->solve("t2_eqns[oo][ff]{u} += 1/2 t2_1[oo][ff]{u} 2@2 <[ff]:[ff]>");

  blas->solve("t2_eqns[oO][fF]{u} +=  t2_1[oO][vV]{u} 2@2 <[ff]|[vv]>");
  blas->solve("t2_eqns[oO][fF]{u} +=  t2_1[oO][vF]{u} 2@2 <[ff]|[vf]>");
  blas->solve("t2_eqns[oO][fF]{u} +=  t2_1[oO][fV]{u} 2@2 <[ff]|[fv]>");
  blas->solve("t2_eqns[oO][fF]{u} +=  t2_1[oO][fF]{u} 2@2 <[ff]|[ff]>");

  blas->solve("t2_eqns[OO][FF]{u} += 1/2 t2_1[OO][VV]{u} 2@2 <[ff]:[vv]>");
  blas->solve("t2_eqns[OO][FF]{u} += 1/2 t2_1[OO][VF]{u} 2@2 <[ff]:[vf]>");
  blas->solve("t2_eqns[OO][FF]{u} += 1/2 t2_1[OO][FV]{u} 2@2 <[ff]:[fv]>");
  blas->solve("t2_eqns[OO][FF]{u} += 1/2 t2_1[OO][FF]{u} 2@2 <[ff]:[ff]>");

  // HPHP Ring
  blas->solve("t2_eqns[oo][ff]{u} += #1342#   t2_1[of][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #1432# - t2_1[of][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #2341# - t2_1[of][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #2431#   t2_1[of][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #1342#   t2_1[of][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #1432# - t2_1[of][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #2341# - t2_1[of][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][ff]{u} += #2431#   t2_1[of][OV]{u} 2@2 ([fo]|[ov])");

  blas->solve("t2_eqns[oo][ff]{u} += #1342#   t2_1[of][of]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #1432# - t2_1[of][of]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #2341# - t2_1[of][of]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #2431#   t2_1[of][of]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #1342#   t2_1[of][OF]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #1432# - t2_1[of][OF]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #2341# - t2_1[of][OF]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oo][ff]{u} += #2431#   t2_1[of][OF]{u} 2@2 ([fo]|[of])");

  blas->solve("t2_eqns[oO][fF]{u} += #1342#   t2_1[of][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oO][fF]{u} += #1342#   t2_1[of][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oO][fF]{u} += #1423# - t2_1[oF][Ov]{u} 2@2 <[of]|[ov]>");
  blas->solve("t2_eqns[oO][fF]{u} += #2314# - t2_1[oV][Of]{u} 1@2 <[of]|[ov]>");
  blas->solve("t2_eqns[oO][fF]{u} += #2431#   t2_1[OF][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oO][fF]{u} += #2431#   t2_1[OF][ov]{u} 2@2 ([fo]:[ov])");

  blas->solve("t2_eqns[oO][fF]{u} += #1342#   t2_1[of][of]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oO][fF]{u} += #1342#   t2_1[of][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[oO][fF]{u} += #1423# - t2_1[oF][Of]{u} 2@2 <[of]|[of]>");
  blas->solve("t2_eqns[oO][fF]{u} += #2314# - t2_1[oF][Of]{u} 1@2 <[of]|[of]>");
  blas->solve("t2_eqns[oO][fF]{u} += #2431#   t2_1[OF][OF]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[oO][fF]{u} += #2431#   t2_1[OF][of]{u} 2@2 ([fo]:[of])");

  blas->solve("t2_eqns[OO][FF]{u} += #1342#   t2_1[OF][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #1432# - t2_1[OF][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #2341# - t2_1[OF][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #2431#   t2_1[OF][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #1342#   t2_1[OF][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #1432# - t2_1[OF][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #2341# - t2_1[OF][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][FF]{u} += #2431#   t2_1[OF][ov]{u} 2@2 ([fo]|[ov])");

  blas->solve("t2_eqns[OO][FF]{u} += #1342#   t2_1[OF][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #1432# - t2_1[OF][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #2341# - t2_1[OF][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #2431#   t2_1[OF][OF]{u} 2@2 ([fo]:[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #1342#   t2_1[OF][of]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #1432# - t2_1[OF][of]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #2341# - t2_1[OF][of]{u} 2@2 ([fo]|[of])");
  blas->solve("t2_eqns[OO][FF]{u} += #2431#   t2_1[OF][of]{u} 2@2 ([fo]|[of])");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][ff]{u}   = t2_eqns[oo][ff]{u} / d2[oo][ff]{u}");
  blas->solve("t2_2[oO][fF]{u}   = t2_eqns[oO][fF]{u} / d2[oO][fF]{u}");
  blas->solve("t2_2[OO][FF]{u}   = t2_eqns[OO][FF]{u} / d2[OO][FF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_2[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff_3 = blas->get_scalar("ECCSD",0);

  fprintf(outfile,"\n      * MP3 (vv)             = %20.12f",E_vv_3);
  fprintf(outfile,"\n      * MP3 (vf)             = %20.12f",E_vf_3);
  fprintf(outfile,"\n      * MP3 (ff)             = %20.12f",E_ff_3);
  fprintf(outfile,"\n      * MP3 total correction = %20.12f",E_vv_3 + E_vf_3 * 2.0 + E_ff_3);

  fflush(outfile);
}








void CCMRCC::perturbative_cbs2()
{

  fprintf(outfile,"\n\n  Computing perturbative corrections for frozen-virtual MOs");

  blas->solve("t2_1[oo][vv]{u}   = <[oo]:[vv]> / d2[oo][vv]{u}");
  blas->solve("t2_1[oO][vV]{u}   = <[oo]|[vv]> / d2[oO][vV]{u}");
  blas->solve("t2_1[OO][VV]{u}   = <[oo]:[vv]> / d2[OO][VV]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     t2_1[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VV]{u} . <[oo]:[vv]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vv = blas->get_scalar("ECCSD",0);

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_1[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  blas->solve("t2_1[oo][vf]{u}   = <[oo]:[vf]> / d2[oo][vf]{u}");
  blas->solve("t2_1[oO][vF]{u}   = <[oo]|[vf]> / d2[oO][vF]{u}");
  blas->solve("t2_1[OO][VF]{u}   = <[oo]:[vf]> / d2[OO][VF]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_1[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf = blas->get_scalar("ECCSD",0);

  blas->solve("t2_1[oo][fv]{u}   = <[oo]:[fv]> / d2[oo][fv]{u}");
  blas->solve("t2_1[oO][fV]{u}   = <[oo]|[fv]> / d2[oO][fV]{u}");
  blas->solve("t2_1[OO][FV]{u}   = <[oo]:[fv]> / d2[OO][FV]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][fv]{u} . <[oo]:[fv]>");
  blas->solve("Eabab{u} =     t2_1[oO][fV]{u} . <[oo]|[fv]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][FV]{u} . <[oo]:[fv]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_fv = blas->get_scalar("ECCSD",0);


  blas->solve("t2_1[oo][ff]{u}   = <[oo]:[ff]> / d2[oo][ff]{u}");
  blas->solve("t2_1[oO][fF]{u}   = <[oo]|[ff]> / d2[oO][fF]{u}");
  blas->solve("t2_1[OO][FF]{u}   = <[oo]:[ff]> / d2[OO][FF]{u}");

  blas->solve("Eaaaa{u} = 1/4 t2_1[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_1[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_1[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff = blas->get_scalar("ECCSD",0);

  fprintf(outfile,"\n\n        CBS second-order correction (vf) = %20.12f",E_vf);
  fprintf(outfile,"\n        CBS second-order correction (fv) = %20.12f",E_fv);
  fprintf(outfile,"\n        CBS second-order correction (ff) = %20.12f",E_ff);
  fprintf(outfile,"\n        CBS second-order correction      = %20.12f",E_ff + E_vf + E_fv);

  // Sort amplitudes
  blas->solve("t2[ov][ov]{u} = #1324# t2[oo][vv]{u}");
  blas->solve("t2[ov][OV]{u} = #1324# t2[oO][vV]{u}");
  blas->solve("t2[oV][Ov]{u} = #1342# t2[oO][vV]{u}");
  blas->solve("t2[OV][OV]{u} = #1324# t2[OO][VV]{u}");

  // OOVF Terms

  // PPPP Ladder
  blas->solve("t2_eqns[oo][vf]{u} += 1/2 t2[oo][vv]{u} 2@2 <[vf]:[vv]>");

  blas->solve("t2_eqns[oO][vF]{u} +=  t2[oO][vV]{u} 2@2 <[vf]|[vv]>");

  blas->solve("t2_eqns[OO][VF]{u} += 1/2 t2[OO][VV]{u} 2@2 <[vf]:[vv]>");

  // HPHP Ring
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2[ov][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2[ov][ov]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #1342#   t2[ov][OV]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oo][vf]{u} += #2341# - t2[ov][OV]{u} 2@2 ([fo]|[ov])");

  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2[ov][ov]{u} 2@2 ([fo]|[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #1342#   t2[ov][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[oO][vF]{u} += #2314# - t2[oV][Ov]{u} 1@2 <[of]|[ov]>");

  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2[OV][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2[OV][OV]{u} 2@2 ([fo]:[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #1342#   t2[ov][OV]{u} 1@2 ([fo]|[ov])");
  blas->solve("t2_eqns[OO][VF]{u} += #2341# - t2[ov][OV]{u} 1@2 ([fo]|[ov])");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][vf]{u}   = t2_eqns[oo][vf]{u} / d2[oo][vf]{u}");
  blas->solve("t2_2[oO][vF]{u}   = t2_eqns[oO][vF]{u} / d2[oO][vF]{u}");
  blas->solve("t2_2[OO][VF]{u}   = t2_eqns[OO][VF]{u} / d2[OO][VF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][vf]{u} . <[oo]:[vf]>");
  blas->solve("Eabab{u} =     t2_2[oO][vF]{u} . <[oo]|[vf]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][VF]{u} . <[oo]:[vf]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_vf_3 = blas->get_scalar("ECCSD",0);


  // OOFF Terms
  // PPPP Ladder
  blas->solve("t2_eqns[oo][ff]{u} += 1/2 t2[oo][vv]{u} 2@2 <[ff]:[vv]>");

  blas->solve("t2_eqns[oO][fF]{u} +=  t2[oO][vV]{u} 2@2 <[ff]|[vv]>");

  blas->solve("t2_eqns[OO][FF]{u} += 1/2 t2[OO][VV]{u} 2@2 <[ff]:[vv]>");

  // Solve for second-order amplitudes
  blas->solve("t2_2[oo][ff]{u}   = t2_eqns[oo][ff]{u} / d2[oo][ff]{u}");
  blas->solve("t2_2[oO][fF]{u}   = t2_eqns[oO][fF]{u} / d2[oO][fF]{u}");
  blas->solve("t2_2[OO][FF]{u}   = t2_eqns[OO][FF]{u} / d2[OO][FF]{u}");

  // Compute third-order energy
  blas->solve("Eaaaa{u} = 1/4 t2_2[oo][ff]{u} . <[oo]:[ff]>");
  blas->solve("Eabab{u} =     t2_2[oO][fF]{u} . <[oo]|[ff]>");
  blas->solve("Ebbbb{u} = 1/4 t2_2[OO][FF]{u} . <[oo]:[ff]>");

  blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");

  double E_ff_3 = blas->get_scalar("ECCSD",0);

  fprintf(outfile,"\n\n        CBS third-order  correction (vf) = %20.12f",E_vf_3);
  fprintf(outfile,"\n        CBS third-order  correction (fv) = %20.12f",E_vf_3);
  fprintf(outfile,"\n        CBS third-order  correction (ff) = %20.12f",E_ff_3);
  fprintf(outfile,"\n        CBS third-order  correction      = %20.12f",E_vf_3 * 2.0 + E_ff_3);

  fprintf(outfile,"\n\n      * CBS corrected energy (2)         = %20.12f",current_energy + E_ff + E_vf + E_fv);
  fprintf(outfile,"\n      * CBS corrected energy (2+3)       = %20.12f",current_energy + E_ff + E_vf + E_fv + E_vf_3 * 2.0 + E_ff_3);
  fflush(outfile);
}

}}

/* End Namespaces */



//// OOFV Terms
//
//blas->solve("t2_eqns[oo][fv]{u}  = 1/2 <[oo]:[oo]> 1@1 t2_1[oo][fv]{u}");
//blas->solve("t2_eqns[oO][fV]{u}  = 1/2 <[oo]|[oo]> 1@1 t2_1[oO][fV]{u}");
//
//blas->solve("t2_eqns[oo][fv]{u} += 1/2 t2_1[oo][vv]{u} 2@2 <[fv]:[vv]>");
//blas->solve("t2_eqns[oo][fv]{u} += 1/2 t2_1[oo][vf]{u} 2@2 <[fv]:[vf]>");
//blas->solve("t2_eqns[oo][fv]{u} += 1/2 t2_1[oo][fv]{u} 2@2 <[fv]:[fv]>");
//blas->solve("t2_eqns[oo][fv]{u} += 1/2 t2_1[oo][ff]{u} 2@2 <[fv]:[ff]>");
//
//blas->solve("t2_eqns[oO][fV]{u} +=  t2_1[oO][vV]{u} 2@2 <[fv]|[vv]>");
//blas->solve("t2_eqns[oO][fV]{u} +=  t2_1[oO][vF]{u} 2@2 <[fv]|[vf]>");
//blas->solve("t2_eqns[oO][fV]{u} +=  t2_1[oO][fV]{u} 2@2 <[fv]|[fv]>");
//blas->solve("t2_eqns[oO][fV]{u} +=  t2_1[oO][fF]{u} 2@2 <[fv]|[ff]>");
//
//blas->solve("t2_2[oo][fv]{u}   = t2_eqns[oo][fv]{u} / d2[oo][fv]{u}");
//blas->solve("t2_2[oO][fV]{u}   = t2_eqns[oO][fV]{u} / d2[oO][fV]{u}");
//blas->solve("t2_2[OO][FV]{u}   = t2_eqns[OO][FV]{u} / d2[OO][FV]{u}");
//
//blas->solve("Eaaaa{u} = 1/4 t2_2[oo][fv]{u} . <[oo]:[fv]>");
//blas->solve("Eabab{u} =     t2_2[oO][fV]{u} . <[oo]|[fv]>");
//blas->solve("Ebbbb{u} = 1/4 t2_2[OO][FV]{u} . <[oo]:[fv]>");
//
//blas->solve("ECCSD{u}  = Eaaaa{u} + Eabab{u} + Ebbbb{u}");
//
//double E_fv_3 = blas->get_scalar("ECCSD",0);
