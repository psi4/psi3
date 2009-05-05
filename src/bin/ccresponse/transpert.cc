/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/* transpert(): Transform various one-electron property integrals from
** the AO to the MO basis.  In some cases, we must also add
** appropriate prefactors and signs.  The only argument is a
** character string indicating the type of integral we want: "Mu",
** "L", "L*", "P", or "P*".  The cints code produces only lower
** triangles, so we must unpack the integrals and keep up with
** symmetric vs. antisymmetric cases.
**
** Notes on specific integrals (all produced by "cints --oeprop") used
** here:
**
** (1) Mu: Length-gauge electric dipole moment integrals = -r.  These
** already include the electronic charge, and they are symmetric wrt
** index permutation.

** (2) L: Magnetic dipole integrals = -i/2 (r x p).  The
** input integrals are angular momentum integrals (r x p) = - i (r x Del), 
** but we multiply these by -0.5 to account for both the sign of the 
** electronic charge and the definition of the magnetic dipole.  These 
** integrals are antisymmetric wrt index permutation.  Use "L*" to select 
** the complex conjugate of the operator (i.e., multiply by -1).
**
** (3) P: Velocity-gauge electric dipole moment integrals = -p.  The input
** integrals are +del integrals, but the definition of
** the linear momentum operator involves a -1 and the electron charge is
**  -1, so no special multiplcation is necessary (i.e., p = -i del, so 
** e*p = i del).  These integrals are antisymmetric wrt to index permutation. 
** Use "P*" to select the complex conjugate of the operator (i.e., multiply 
** by -1).
**
** NB: The magnetic-dipole and velocity-gauge electric-dipole operators are
** both pure-imaginary, which means that one must take care to account for
** the factor of i included implicity in their definition.  This matters,
** for example, in the computation of optical rotations.  See the notes in
** optrot.cc for specifics.
**
** -TDC, 11/05
** Updated by TDC 4/09
*/

void transpert(const char *pert)
{
  int nao, nso, nmo, noei_ao;
  int alpha;
  int i, j, ij;
  double *scratch, **TMP, **X, **target;
  const char *name;
  double prefactor, anti, sign;

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = moinfo.noei_ao;

  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);

  if(!strcmp(pert,"Mu")) { prefactor = 1.0; anti = 1.0; sign = 1.0; }
  else if(!strcmp(pert, "L")) { prefactor = -0.5; anti = -1.0; sign = 1.0; }
  else if(!strcmp(pert, "L*")) { prefactor = -0.5; anti = -1.0; sign = -1.0; }
  else if(!strcmp(pert, "P")) { prefactor = 1.0; anti = -1.0; sign = 1.0; }
  else if(!strcmp(pert, "P*")) { prefactor = 1.0; anti = -1.0; sign = -1.0; }

  for(alpha=0; alpha < 3; alpha++) {

    target = block_matrix(nmo,nmo);

    if(!strcmp(pert,"Mu")) {
      if(alpha == 0) { name = PSIF_AO_MX; moinfo.MUX = target; }
      else if(alpha == 1) { name = PSIF_AO_MY; moinfo.MUY = target; }
      else if(alpha == 2) { name = PSIF_AO_MZ; moinfo.MUZ = target; }
    }
    else if(!strcmp(pert,"L") || !strcmp(pert, "L*")) {
      if(alpha == 0) { name = PSIF_AO_LX; moinfo.LX = target; }
      else if(alpha == 1) { name = PSIF_AO_LY; moinfo.LY = target; }
      else if(alpha == 2) { name = PSIF_AO_LZ; moinfo.LZ = target; }
    }
    else if(!strcmp(pert,"P") || !strcmp(pert, "P*")) {
      if(alpha == 0) { name = PSIF_AO_NablaX; moinfo.PX = target; }
      else if(alpha == 1) { name = PSIF_AO_NablaY; moinfo.PY = target; }
      else if(alpha == 2) { name = PSIF_AO_NablaZ; moinfo.PZ = target; }
    }

    iwl_rdone(PSIF_OEI, name, scratch, noei_ao, 0, 0, outfile);
    for(i=0,ij=0; i < nao; i++)
      for(j=0; j <= i; j++,ij++) {
	TMP[i][j] = prefactor * sign * scratch[ij];
	TMP[j][i] = anti * prefactor * sign * scratch[ij];
      }

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
	    0,&(X[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
	    0,&(TMP[0][0]),nao);

    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
	    0,&(X[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
	    0,&(target[0][0]),nmo);

    zero_arr(scratch,noei_ao);

  }

  free(scratch);
  free_block(TMP);
  free_block(X);
}

}} // namespace psi::ccresponse
