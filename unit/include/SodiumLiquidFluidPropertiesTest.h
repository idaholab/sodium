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
