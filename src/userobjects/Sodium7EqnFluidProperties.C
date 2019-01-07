#include "Sodium7EqnFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

void DIFF_ps_t_Na(double t, double & ps, double & dpsdt, double & d2psdt2);
int DIFF_ts_p_Na(double p, double & ts, double & dtsdp, double & d2tsdp2);

const Real Sodium7EqnFluidProperties::_P_critical = 25.64E+6;

registerMooseObject("SodiumApp", Sodium7EqnFluidProperties);

template <>
InputParameters
validParams<Sodium7EqnFluidProperties>()
{
  InputParameters params = validParams<TwoPhaseFluidProperties>();
  params += validParams<NaNInterface>();
  params.addClassDescription("Fluid properties of sodium for the 7-equation model.");
  return params;
}

Sodium7EqnFluidProperties::Sodium7EqnFluidProperties(const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),
    NaNInterface(this),
    _to_atm(1. / 101325.),
    _to_Pa(101325.),
    _to_R(9. / 5.),
    _to_K(5. / 9.)
{
  {
    std::string class_name = "SodiumLiquidFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, _liquid_name, params);
  }
  _fp_liquid = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_liquid_name);

  {
    std::string class_name = "SodiumVaporFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    _fe_problem.addUserObject(class_name, _vapor_name, params);
  }
  _fp_vapor = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_vapor_name);
}

Real
Sodium7EqnFluidProperties::p_critical() const
{
  return _P_critical;
}

Real
Sodium7EqnFluidProperties::T_sat(Real pressure) const
{
  pressure *= _to_atm;

  static const double p0 = 1.E-8 * 0.101325; // just a lower limit, not actual triple point press.
  static const double pc = 25.64 * 0.101325;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    DIFF_ts_p_Na(pressure, ts, dtsdp, d2tsdp2);
    return ts * _to_K;
  }
  else
    return getNaN();
}

Real
Sodium7EqnFluidProperties::p_sat(Real temperature) const
{
  temperature *= _to_R;

  static const double t0 = 273.15 * (9. / 5.); // just a lower limit, not actual triple point temp.
  static const double tc = 2503.7 * (9. / 5.);

  if (t0 < temperature && temperature < tc)
  {
    double ps, dpsdt, d2psdt2;
    DIFF_ps_t_Na(temperature, ps, dpsdt, d2psdt2);
    return ps * _to_Pa;
  }
  else
    return getNaN();
}

Real
Sodium7EqnFluidProperties::dT_sat_dp(Real pressure) const
{
  pressure *= _to_atm;

  static const double p0 = 1.E-8 * 0.101325; // just a lower limit, not actual triple point press.
  static const double pc = 25.64 * 0.101325;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    DIFF_ts_p_Na(pressure, ts, dtsdp, d2tsdp2);
    return dtsdp * _to_K / _to_Pa;
  }
  else
    return getNaN();
}
