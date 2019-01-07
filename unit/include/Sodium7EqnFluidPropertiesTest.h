#ifndef SODIUM7EQNFLUIDPROPERTIESTEST_H
#define SODIUM7EQNFLUIDPROPERTIESTEST_H

#include "MooseObjectUnitTest.h"
#include "Sodium7EqnFluidProperties.h"

class Sodium7EqnFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  Sodium7EqnFluidPropertiesTest() : MooseObjectUnitTest("SodiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("Sodium7EqnFluidProperties");
    _fe_problem->addUserObject("Sodium7EqnFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<Sodium7EqnFluidProperties>("fp");
  }

  const Sodium7EqnFluidProperties * _fp;
};

#endif /* SODIUM7EQNFLUIDPROPERTIESTEST_H */
