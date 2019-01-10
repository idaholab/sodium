#include "SodiumLiquidFluidPropertiesTest.h"
#include "SinglePhaseFluidPropertiesTestUtils.h"

TEST_F(SodiumLiquidFluidPropertiesTest, test)
{
  const Real T = 1500.;
  const Real p = 1101124.69860416;

  const Real rho_from_p_T = _fp->rho_from_p_T(p, T);
  const Real rho = rho_from_p_T;

  const Real h_from_p_T = _fp->h_from_p_T(p, T);
  const Real h = h_from_p_T;

  const Real e_from_p_rho = _fp->e_from_p_rho(p, rho);
  const Real e = e_from_p_rho;

  const Real v = 1 / rho;

  const Real s_from_v_e = _fp->s_from_v_e(v, e);
  const Real s = s_from_v_e;

  // p
  REL_TEST(_fp->p_from_v_e(v, e), p, REL_TOL_CONSISTENCY);
  REL_TEST(_fp->p_from_h_s(h, s), p, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->p_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->p_from_h_s, h, s, REL_TOL_DERIVATIVE);

  // T
  REL_TEST(_fp->T_from_v_e(v, e), T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->T_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // rho
  REL_TEST(rho_from_p_T, 656.62644694597589, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->rho_from_p_s(p, s), rho_from_p_T, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->rho_from_p_T, p, T, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->rho_from_p_s, p, s, REL_TOL_DERIVATIVE);

  // e
  REL_TEST(e_from_p_rho, 1959009.7928613776, REL_TOL_SAVED_VALUE);
  REL_TEST(_fp->e_from_v_h(v, h), e, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->e_from_p_rho, p, rho, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->e_from_v_h, v, h, REL_TOL_DERIVATIVE);

  // c
  const Real c = _fp->c_from_v_e(v, e);
  REL_TEST(c, 1758.5592982054891, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->c_from_v_e, v, e, REL_TOL_DERIVATIVE);

  // cp
  const Real cp = _fp->cp_from_v_e(v, e);
  REL_TEST(cp, 1417.2885128554665, REL_TOL_SAVED_VALUE);

  // cv
  const Real cv = _fp->cv_from_v_e(v, e);
  REL_TEST(cv, 978.54515799090188, REL_TOL_SAVED_VALUE);

  // mu
  const Real mu = _fp->mu_from_v_e(v, e);
  REL_TEST(mu, 0.00012410852760561553, REL_TOL_SAVED_VALUE);

  // k
  const Real k = _fp->k_from_v_e(v, e);
  REL_TEST(k, 39.339228017959506, REL_TOL_SAVED_VALUE);

  // s
  REL_TEST(s, 4646.959555608144, REL_TOL_EXTERNAL_VALUE);
  REL_TEST(_fp->s_from_h_p(h, p), s, REL_TOL_CONSISTENCY);
  DERIV_TEST(_fp->s_from_v_e, v, e, REL_TOL_DERIVATIVE);
  DERIV_TEST(_fp->s_from_h_p, h, p, REL_TOL_DERIVATIVE);

  // g
  REL_TEST(_fp->g_from_v_e(v, e), -5009752.5981984176, REL_TOL_EXTERNAL_VALUE);

  // h
  REL_TEST(h_from_p_T, 1960686.7352137996, REL_TOL_SAVED_VALUE);
  DERIV_TEST(_fp->h_from_p_T, p, T, REL_TOL_DERIVATIVE);

  // sigma
  REL_TEST(_fp->sigma_from_p_T(p, T), 0.085924584655757016, REL_TOL_SAVED_VALUE);
}
