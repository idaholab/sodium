#include "Sodium7EqnFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(Sodium7EqnFluidPropertiesTest, test)
{
  Real T = 1200.;  // K
  Real p = 101325; // Pa

  // Tsat + derivatives
  REL_TEST(_fp->T_sat(p), 1154.6031110050355, REL_TOL_SAVED_VALUE);
  DERIV_TEST_1D(_fp->T_sat, _fp->dT_sat_dp, p, REL_TOL_DERIVATIVE);

  // Psat
  REL_TEST(_fp->p_sat(T), 150357.12693001426, REL_TOL_SAVED_VALUE);

  // Psat
  REL_TEST(_fp->sigma_from_T(T), 0.11534571734609324, REL_TOL_SAVED_VALUE);
  DERIV_TEST_1D(_fp->sigma_from_T, _fp->dsigma_dT_from_T, T, REL_TOL_DERIVATIVE);
}
