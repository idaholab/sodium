#pragma once

#include "MooseApp.h"

class Factory;
class SodiumApp;

template <>
InputParameters validParams<SodiumApp>();

class SodiumApp : public MooseApp
{
public:
  SodiumApp(InputParameters parameters);

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjects(Factory & factory);

protected:
};
