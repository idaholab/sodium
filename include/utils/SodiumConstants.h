//* This file is part of sodium
//* https://github.com/idaholab/sodium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/sodium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

namespace SodiumConstants
{
  // Note that the molar mass can change based on the distribution of monomers,
  // dimers, tetramers, etc. For now, the number used here is the molar weight
  // one gets from a periodic table; it is the molar mass associated with only
  // monomers, with a natural isotopic distribution.
  constexpr Real molar_mass = 22.989769e-3;

  // Values taken from Sodium-NaK Engineering Handbook, Volume 1
  constexpr Real critical_temperature = 2733.0;
  constexpr Real critical_pressure = 41.3406e6;
}
