#include "SodiumTwoPhaseFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(SodiumTwoPhaseFluidPropertiesTest, test)
{
  Real T = 1200.;  // K
  Real p = 101325; // Pa

  // T_triple
  REL_TEST(_fp->T_triple(), 370.98, REL_TOL_SAVED_VALUE);

  // L_fusion
  REL_TEST(_fp->L_fusion(), 0.112968, REL_TOL_SAVED_VALUE);

  // Tsat
  REL_TEST(_fp->T_sat(_fp->p_sat(371.0)), 371.0, REL_TOL_CONSISTENCY);
  REL_TEST(_fp->T_sat(p), 1154.6031110050355, REL_TOL_SAVED_VALUE);
  DERIV_TEST_1D(_fp->T_sat, _fp->dT_sat_dp, p, REL_TOL_DERIVATIVE);

  // Psat
  REL_TEST(_fp->p_sat(T), 150357.12693001426, REL_TOL_SAVED_VALUE);

  // sigma
  REL_TEST(_fp->sigma_from_T(T), 0.11534571734609324, REL_TOL_SAVED_VALUE);
  DERIV_TEST_1D(_fp->sigma_from_T, _fp->dsigma_dT_from_T, T, REL_TOL_DERIVATIVE);
}
