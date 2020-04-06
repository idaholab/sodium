#pragma once

#include "MooseApp.h"

class Factory;

class SodiumApp : public MooseApp
{
public:
  SodiumApp(InputParameters parameters);

public:
  static InputParameters validParams();
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjects(Factory & factory);
};
