//* This file is part of sodium
//* https://github.com/idaholab/sodium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/sodium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SodiumApp.h"
#include "SodiumRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "FluidPropertiesApp.h"

InputParameters
SodiumApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

registerKnownLabel("SodiumApp");

SodiumApp::SodiumApp(InputParameters parameters)
  : MooseApp(parameters)
{
  SodiumApp::registerAll(_factory, _action_factory, _syntax);
}

// External entry point for dynamic application loading
extern "C" void
SodiumApp__registerApps()
{
  SodiumApp::registerApps();
}

void
SodiumApp::registerApps()
{
  registerApp(SodiumApp);
}

extern "C" void
SodiumApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SodiumApp::registerAll(f, af, s);
}

void
SodiumApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"SodiumApp"});
  Registry::registerActionsTo(af, {"SodiumApp"});

  FluidPropertiesApp::registerAll(f, af, s);
}

void
SodiumApp::registerObjects(Factory & factory)
{
  mooseDeprecated("SodiumApp: use registerAll instead of registerObjects");
  Registry::registerObjectsTo(factory, {"SodiumApp"});
}
