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
