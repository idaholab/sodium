#ifndef SODIUM7EQNFLUIDPROPERTIES_H
#define SODIUM7EQNFLUIDPROPERTIES_H

#include "TwoPhaseFluidProperties.h"
#include "NaNInterface.h"

class Sodium7EqnFluidProperties;
class SinglePhaseFluidProperties;

template <>
InputParameters validParams<Sodium7EqnFluidProperties>();

/**
 * Sodium interface for 7-eqn model
 */
class Sodium7EqnFluidProperties : public TwoPhaseFluidProperties, public NaNInterface
{
public:
  Sodium7EqnFluidProperties(const InputParameters & parameters);

  virtual Real p_critical() const;
  virtual Real T_sat(Real pressure) const;
  virtual Real p_sat(Real temperature) const;
  virtual Real dT_sat_dp(Real pressure) const;

protected:
  // Critical pressure
  static const Real _P_critical;

protected:
  /// Conversion factor from Pa to atm
  const Real _to_atm;
  /// Conversion factor from atm to Pa
  const Real _to_Pa;
  /// Conversion factor from K to R
  const Real _to_R;
  /// Conversion factor from R to K
  const Real _to_K;
};

#endif /* SODIUM7EQNFLUIDPROPERTIES_H */
