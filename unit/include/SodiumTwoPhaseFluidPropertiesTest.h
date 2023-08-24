//* This file is part of sodium
//* https://github.com/idaholab/sodium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/sodium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObjectUnitTest.h"
#include "SodiumTwoPhaseFluidProperties.h"

class SodiumTwoPhaseFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  SodiumTwoPhaseFluidPropertiesTest() : MooseObjectUnitTest("SodiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("SodiumTwoPhaseFluidProperties");
    _fe_problem->addUserObject("SodiumTwoPhaseFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<SodiumTwoPhaseFluidProperties>("fp");
  }

  const SodiumTwoPhaseFluidProperties * _fp;
};
