/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include <libutil/libutil.h>

#include "blas.h"
#include "debugging.h"
#include "idmrpt2.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void IDMRPT2::build_F_intermediates()
{
  build_F_ae_intermediates();
  build_F_AE_intermediates();

  build_F_mi_intermediates();
  build_F_MI_intermediates();
}

/**
 * \brief Computes the intermediate
 * \f[
 * \mathcal{F}_{ae}(\mu) = (1-\delta_{ae}) < a| \hat{F}_{\mu}^{A} + \hat{F}_{\mu}^{E}|e>
 * \f]
 */
void IDMRPT2::build_F_ae_intermediates()
{
  START_TIMER(1,"Building the F_ae Intermediates");

  blas->solve("F_ae[v][v]{u} = fock[v][v]{u}");
  blas->solve_zero_two_diagonal("F_ae[v][v]{u}");
  blas->zero_non_external("F_ae[v][v]{u}");

  DEBUGGING(3,
    blas->print("F_ae[v][v]{u}");
  );
  END_TIMER(1);
}

/**
 * \brief Computes the intermediate
 * \f[
 * \mathcal{F}_{AE}(\mu) = (1-\delta_{AE}) < A| \hat{F}_{\mu}^{A} + \hat{F}_{\mu}^{E}|E>
 * \f]
 */
void IDMRPT2::build_F_AE_intermediates()
{
  START_TIMER(1,"Building the F_AE Intermediates");

  blas->solve("F_AE[V][V]{u} = fock[V][V]{u}");
  blas->solve_zero_two_diagonal("F_AE[V][V]{u}");
  blas->zero_non_external("F_AE[V][V]{u}");

  DEBUGGING(3, blas->print("F_AE[V][V]{u}"); );
  END_TIMER(1);
}

/**
 * \brief Computes the intermediate
 * \f[
 * \mathcal{F}_{mi}(\mu) = (1-\delta_{mi}) <m| \hat{F}_{\mu}^{D} + \hat{F}_{\mu}^{A}|i>
 * \f]
 */
void IDMRPT2::build_F_mi_intermediates()
{
  START_TIMER(1,"Building the F_mi Intermediates");

  blas->solve("F_mi[o][o]{u} = fock[o][o]{u}");
  blas->solve_zero_two_diagonal("F_mi[o][o]{u}");
  blas->zero_non_doubly_occupied("F_mi[o][o]{u}");

  DEBUGGING(3, blas->print("F_mi[o][o]{u}"); );
  END_TIMER(1);
}

/**
 * \brief Computes the intermediate
 * \f[
 * \mathcal{F}_{MI}(\mu) = (1-\delta_{MI}) <M| \hat{F}_{\mu}^{D} + \hat{F}_{\mu}^{A}|I>
 * \f]
 */
void IDMRPT2::build_F_MI_intermediates()
{
  START_TIMER(1,"Building the F_MI Intermediates");

  blas->solve("F_MI[O][O]{u} = fock[O][O]{u}");
  blas->solve_zero_two_diagonal("F_MI[O][O]{u}");
  blas->zero_non_doubly_occupied("F_MI[O][O]{u}");

  DEBUGGING(3,
    blas->print("F_MI[O][O]{u}");
  );
  END_TIMER(1);
}

}} /* End Namespaces */
