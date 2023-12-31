//* This file is part of sodium
//* https://github.com/idaholab/sodium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/sodium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SodiumTwoPhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"
#include "contrib/libSodiumProperties/Na_Golden.h"

const Real SodiumTwoPhaseFluidProperties::_P_critical = 25.64E+6;

// Value is taken from NIST Chemistry WebBook, SRD 69.
// Original reference:
// R. E. Honing and D. A. Kramer. Vapor pressure data for the solid and liquid elements (1969)
const Real SodiumTwoPhaseFluidProperties::_T_triple = 370.98;

// Value is taken from the following reference:
//
// O. J. Foust. Sodium-NaK Engineering Handbook, Volume 1: Sodium Chemistry and Physical Properties
// (1972). Division of Reactor Development and Technology, United States Atomic Energy Commission.
//
// Value given was 27 cal/g, which was converted to J/kg.
const Real SodiumTwoPhaseFluidProperties::_L_fusion = 112968.0;

registerMooseObject("SodiumApp", SodiumTwoPhaseFluidProperties);

InputParameters
SodiumTwoPhaseFluidProperties::validParams()
{
  InputParameters params = TwoPhaseFluidProperties::validParams();
  params += NaNInterface::validParams();
  params.addClassDescription("Two-phase sodium fluid properties");
  return params;
}

SodiumTwoPhaseFluidProperties::SodiumTwoPhaseFluidProperties(const InputParameters & parameters)
  : TwoPhaseFluidProperties(parameters),
    NaNInterface(this),
    _to_atm(1. / 101325.),
    _to_Pa(101325.),
    _to_R(9. / 5.),
    _to_K(5. / 9.)
{
  if (_tid == 0)
  {
    std::string class_name = "SodiumLiquidFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    params.set<MooseEnum>("emit_on_nan") = getParam<MooseEnum>("emit_on_nan");
    _fe_problem.addUserObject(class_name, _liquid_name, params);
  }
  _fp_liquid = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_liquid_name, _tid);

  if (_tid == 0)
  {
    std::string class_name = "SodiumVaporFluidProperties";
    InputParameters params = _app.getFactory().getValidParams(class_name);
    params.set<MooseEnum>("emit_on_nan") = getParam<MooseEnum>("emit_on_nan");
    _fe_problem.addUserObject(class_name, _vapor_name, params);
  }
  _fp_vapor = &_fe_problem.getUserObject<SinglePhaseFluidProperties>(_vapor_name, _tid);
}

Real
SodiumTwoPhaseFluidProperties::p_critical() const
{
  return _P_critical;
}

Real
SodiumTwoPhaseFluidProperties::T_triple() const
{
  return _T_triple;
}

Real
SodiumTwoPhaseFluidProperties::L_fusion() const
{
  return _L_fusion;
}

Real
SodiumTwoPhaseFluidProperties::T_sat(Real pressure) const
{
  pressure *= _to_atm;

  // Limits taken as p_sat @ T = 350 K and T = 2500 K
  static const double p0 = 1.5210548366e-06 * _to_atm;
  static const double pc = 2.4232774980e+07 * _to_atm;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    int ierr = DIFF_ts_p_Na(pressure, ts, dtsdp, d2tsdp2);
    if (ierr != 0)
      return getNaN();
    else
      return ts * _to_K;
  }
  else
    return getNaN();
}

Real
SodiumTwoPhaseFluidProperties::p_sat(Real temperature) const
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
SodiumTwoPhaseFluidProperties::dT_sat_dp(Real pressure) const
{
  pressure *= _to_atm;

  // Limits taken as p_sat @ T = 350 K and T = 2500 K
  static const double p0 = 1.5210548366e-06 * _to_atm;
  static const double pc = 2.4232774980e+07 * _to_atm;

  if (p0 <= pressure && pressure <= pc)
  {
    double ts, dtsdp, d2tsdp2;
    int ierr = DIFF_ts_p_Na(pressure, ts, dtsdp, d2tsdp2);
    if (ierr != 0)
      return getNaN();
    else
      return dtsdp * _to_K / _to_Pa;
  }
  else
    return getNaN();
}

Real
SodiumTwoPhaseFluidProperties::sigma_from_T(Real T) const
{
  double sigma, dsigmadt, d2sigmadt2;
  sigma_t_Na(T, sigma, dsigmadt, d2sigmadt2);
  return sigma * 1e-3;
}

Real
SodiumTwoPhaseFluidProperties::dsigma_dT_from_T(Real T) const
{
  double sigma, dsigmadt, d2sigmadt2;
  sigma_t_Na(T, sigma, dsigmadt, d2sigmadt2);
  return dsigmadt * 1e-3;
}
