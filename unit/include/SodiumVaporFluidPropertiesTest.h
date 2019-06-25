#pragma once

#include "MooseObjectUnitTest.h"
#include "SodiumVaporFluidProperties.h"

class SodiumVaporFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  SodiumVaporFluidPropertiesTest() : MooseObjectUnitTest("SodiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("SodiumVaporFluidProperties");
    _fe_problem->addUserObject("SodiumVaporFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObjectTempl<SodiumVaporFluidProperties>("fp");
  }

  const SodiumVaporFluidProperties * _fp;
};
