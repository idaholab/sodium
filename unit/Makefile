###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# HERD_TRUNK_DIR   - Location of the HERD repository
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
SODIUM_DIR        ?= $(shell dirname `pwd`)
MOOSE_SUBMODULE    := $(SODIUM_DIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname $(SODIUM_DIR))/moose
endif
# check that MOOSE is available
MOOSE_CONTENT      := $(shell ls $(MOOSE_DIR) 2> /dev/null)
ifeq ($(MOOSE_CONTENT),)
  $(error MOOSE framework does not seem to be available. Make sure that either the submodule is checked out or that your MOOSE_DIR points to the correct location)
endif

MOOSE_DIR          ?= $(SODIUM_DIR)/moose
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
###############################################################################
CURRENT_DIR        := $(shell pwd)

# Extra stuff for GTEST
ADDITIONAL_INCLUDES += -I$(FRAMEWORK_DIR)/contrib/gtest
ADDITIONAL_LIBS     += $(FRAMEWORK_DIR)/contrib/gtest/libgtest.la

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
################################## MODULES ####################################
FLUID_PROPERTIES  := yes
include           $(MOOSE_DIR)/modules/modules.mk
###############################################################################

include $(SODIUM_DIR)/libSodiumProperties.mk

# dep apps
APPLICATION_DIR    := $(SODIUM_DIR)
APPLICATION_NAME   := sodium
include            $(FRAMEWORK_DIR)/app.mk

APPLICATION_DIR    := $(CURRENT_DIR)
APPLICATION_NAME   := sodium-unit
BUILD_EXEC         := yes
app_BASE_DIR       :=      # Intentionally blank
DEP_APPS           ?= $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include $(FRAMEWORK_DIR)/app.mk

# Find all the sodium unit test source files and include their dependencies.
sodium-unit_srcfiles := $(shell find $(CURRENT_DIR) -name "*.C")
sodium-unit_deps := $(patsubst %.C, %.$(obj-suffix).d, $(sodium-unit_srcfiles))
-include $(sodium-unit_deps)

###############################################################################
# Additional special case targets should be added here
