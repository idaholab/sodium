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
#include "SodiumLiquidFluidProperties.h"

class SodiumLiquidFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  SodiumLiquidFluidPropertiesTest() : MooseObjectUnitTest("SodiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("SodiumLiquidFluidProperties");
    _fe_problem->addUserObject("SodiumLiquidFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<SodiumLiquidFluidProperties>("fp");
  }

  const SodiumLiquidFluidProperties * _fp;
};
