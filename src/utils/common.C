#include "Na_Golden.h"
#include "common.h"

Real
T_sat(Real pressure)
{
  pressure /= 101325.; // to atm

  double ts, dtsdp, d2tsdp2;
  DIFF_ts_p_Na(pressure, ts, dtsdp, d2tsdp2);
  return ts * 5. / 9.; // to K
}

Real
p_sat(Real temperature)
{
  temperature*=9. / 5.; // to R

  double ps, dpsdt, d2psdt2;
  DIFF_ps_t_Na(temperature, ps, dpsdt, d2psdt2);
  return ps*101325.; // to Pa
}

Real
dT_dP_sat(Real pressure)
{
  pressure /= 101325.; // to atm

  double ts, dtsdp, d2tsdp2;
  DIFF_ts_p_Na(pressure, ts, dtsdp, d2tsdp2);
  return dtsdp * 5. / 9. /101325.; // to K/Pa
}

Real
sigma(Real temperature)
{
  double sigma, dsigmadt, d2sigmadt2;
  sigma_t_Na(temperature, sigma, dsigmadt, d2sigmadt2);
  return sigma * 1e-3;
}
