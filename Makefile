# Build O2FoCalGenQC.so for: #pragma cling load("./O2FoCalGenQC.so")
#
# Expects an aliBuild O2 installation. O2_ROOT / ROOTSYS / BOOST_ROOT / VMC_ROOT
# are used when set (alienv); otherwise they are taken from $HOME/alice/sw.

HOME_DIR := $(shell echo $$HOME)

ifeq ($(origin O2_ROOT), undefined)
  O2_ROOT := $(firstword $(wildcard $(HOME_DIR)/alice/sw/*/O2/latest))
endif
ifeq ($(O2_ROOT),)
  $(error O2_ROOT is not set and no $(HOME)/alice/sw/*/O2/latest was found)
endif

ALICE_ARCH := $(abspath $(O2_ROOT)/../..)

ifeq ($(origin ROOTSYS), undefined)
  ROOTSYS := $(ALICE_ARCH)/ROOT/latest
endif
ifeq ($(origin BOOST_ROOT), undefined)
  BOOST_ROOT := $(ALICE_ARCH)/boost/latest
endif
ifeq ($(origin VMC_ROOT), undefined)
  VMC_ROOT := $(ALICE_ARCH)/VMC/latest
endif

ROOT_CONFIG := $(shell command -v root-config 2>/dev/null)
ifeq ($(ROOT_CONFIG),)
  ROOT_CONFIG := $(ROOTSYS)/bin/root-config
endif
ifeq ($(wildcard $(ROOT_CONFIG)),)
  $(error root-config not found. Enter the O2 environment or set ROOTSYS)
endif
ifeq ($(wildcard $(BOOST_ROOT)/include/boost/container_hash/hash.hpp),)
  $(error Boost headers not found under BOOST_ROOT=$(BOOST_ROOT))
endif
ifeq ($(wildcard $(VMC_ROOT)/include/vmc/TMCProcess.h),)
  $(error VMC headers not found under VMC_ROOT=$(VMC_ROOT))
endif
ifeq ($(wildcard $(O2_ROOT)/include/GPU/GPUCommonRtypes.h),)
  $(error O2 headers not found under O2_ROOT=$(O2_ROOT))
endif

# Make's built-in CXX is g++, so ?= cannot override it. Prefer the aliBuild
# compiler that built this O2/ROOT, unless CXX was set on the command line
# or in the environment.
ALICE_CXX := $(ALICE_ARCH)/GCC-Toolchain/latest/bin/c++
ifeq ($(origin CXX), default)
  ifneq ($(wildcard $(ALICE_CXX)),)
    CXX := $(ALICE_CXX)
  else
    CXX := $(shell $(ROOT_CONFIG) --cxx)
  endif
endif

ROOT_CFLAGS := $(shell $(ROOT_CONFIG) --cflags)
ROOT_LIBS   := $(shell $(ROOT_CONFIG) --libs)

PKG_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
TARGET  := $(PKG_DIR)/O2FoCalGenQC.so
SRC     := $(PKG_DIR)/src/O2FoCalGenQC.cpp
OBJ     := $(PKG_DIR)/O2FoCalGenQC.o
DEP     := $(OBJ:.o=.d)

CXXFLAGS := -O2 -fPIC -Wall -Wextra $(ROOT_CFLAGS)
INCLUDES := -I$(PKG_DIR)/src \
            -isystem $(O2_ROOT)/include \
            -isystem $(O2_ROOT)/include/GPU \
            -isystem $(VMC_ROOT)/include/vmc \
            -isystem $(BOOST_ROOT)/include

O2_LIBS := -L$(O2_ROOT)/lib -Wl,-rpath,$(O2_ROOT)/lib \
           -lO2FOCALBase -lO2SimulationDataFormat -lO2CommonUtils

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -shared -fPIC -Wl,-soname,$(notdir $(TARGET)) -o $@ $^ $(O2_LIBS) $(ROOT_LIBS)

$(OBJ): $(SRC) $(PKG_DIR)/src/O2FoCalGenQC.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $(SRC) -o $@

-include $(DEP)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)
