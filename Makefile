# FlatLand — dependency-free build.
#
# No third-party libraries: just a C++17 compiler with threads.
#
#   make              # build ./flatland (the CLI)
#   make lib          # build the static and shared libraries
#   make all          # both
#   make test         # build, then run the whole test suite
#   make install      # install to $(PREFIX), default /usr/local
#   make clean        # remove build artifacts
#
# Override the compiler if needed:  make CXX=clang++

CXX      ?= c++
CC       ?= cc
CXXFLAGS ?= -std=c++17 -O2 -Wall
LDFLAGS  ?= -pthread

# Determinism, not taste, and deliberately NOT part of CXXFLAGS so that
# overriding those cannot silently drop it.
#
# -ffp-contract=off forbids the compiler from fusing a multiply and an add into
# a single FMA. The rasterizer needs that: it decides coverage with the three
# edge functions, and relies on edge(a,b,p) == -edge(b,a,p) holding EXACTLY, so
# a pixel centre lying on the shared edge of two triangles is claimed by both
# rather than by neither. Contraction rounds the two sides differently, the
# antisymmetry breaks, and one-pixel cracks open along shared edges.
#
# This is invisible on baseline x86-64, which has no FMA instruction, and shows
# up on arm64, which does: a subdivided unit cube at -r 0.1 covers 100 pixels on
# one and 98 on the other. A tool whose whole output is a measurement has to give
# the same answer on every machine, so the fusion goes.
# Probed via stdin rather than /dev/null, which MSYS2 maps onto NUL.
FPFLAGS  := $(shell echo | $(CXX) -ffp-contract=off -E -x c++ - >/dev/null 2>&1 \
                    && echo -ffp-contract=off)
PREFIX   ?= /usr/local
PYTHON   ?= python3

INCLUDES  = -Isrc -Iinclude

# The engine. Shared by the CLI, the C API and anything else that links it.
CORE_SRC = src/fl_mesh.cpp \
           src/fl_deflate.cpp \
           src/fl_parallel.cpp \
           src/fl_io.cpp \
           src/fl_field.cpp \
           src/fl_raster.cpp \
           src/fl_project.cpp \
           src/fl_image.cpp \
           src/fl_batch.cpp
CAPI_SRC = src/capi.cpp
CLI_SRC  = cli/main.cpp

CORE_OBJ = $(CORE_SRC:.cpp=.o)
CAPI_OBJ = $(CAPI_SRC:.cpp=.o)
CLI_OBJ  = $(CLI_SRC:.cpp=.o)

# Position-independent copies for the shared library.
LIB_PIC  = $(CORE_SRC:.cpp=.pic.o) $(CAPI_SRC:.cpp=.pic.o)

DEP      = $(CORE_OBJ:.o=.d) $(CAPI_OBJ:.o=.d) $(CLI_OBJ:.o=.d) $(LIB_PIC:.o=.d)

STATICLIB = libflatland.a

# Platform naming. macOS wants .dylib and -install_name; Windows (MinGW/MSYS2/
# Cygwin) wants a .dll plus a separate import library, has no -soname, and puts
# an .exe suffix on executables whether asked to or not — so the CLI target has
# to carry that suffix, otherwise make never sees the file it just built.
UNAME_S := $(shell uname -s)
EXE       =
IMPLIB    =
ifeq ($(UNAME_S),Darwin)
  SHAREDLIB = libflatland.dylib
  SOFLAGS   = -dynamiclib -install_name @rpath/$(SHAREDLIB)
else ifneq (,$(filter MINGW% MSYS% CYGWIN%,$(UNAME_S)))
  EXE       = .exe
  SHAREDLIB = flatland.dll
  IMPLIB    = libflatland.dll.a
  SOFLAGS   = -shared -Wl,--out-implib,$(IMPLIB)
else
  SHAREDLIB = libflatland.so
  SOFLAGS   = -shared -Wl,-soname,$(SHAREDLIB)
endif

BIN       = flatland$(EXE)

.PHONY: all lib test validate install clean

$(BIN): $(CORE_OBJ) $(CLI_OBJ)
	$(CXX) $(CXXFLAGS) $(FPFLAGS) $^ -o $@ $(LDFLAGS)

all: $(BIN) lib
lib: $(STATICLIB) $(SHAREDLIB)

$(STATICLIB): $(CORE_OBJ) $(CAPI_OBJ)
	$(AR) rcs $@ $^

# -fvisibility=hidden keeps the engine's C++ symbols internal, so the shared
# library exports only what flatland.h marks FL_API.
$(SHAREDLIB): $(LIB_PIC)
	$(CXX) $(CXXFLAGS) $(FPFLAGS) $(SOFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(FPFLAGS) $(INCLUDES) $(FL_API_DEF) -MMD -MP -c $< -o $@

%.pic.o: %.cpp
	$(CXX) $(CXXFLAGS) $(FPFLAGS) $(INCLUDES) $(FL_API_DEF) -fPIC -fvisibility=hidden -MMD -MP -c $< -o $@

# The C ABI's only Windows-visible difference between the two libraries is which
# way FL_API points, and that is decided per object file, not per source file.
# The plain .o goes into the static archive (no dllexport); the .pic.o goes into
# the shared library (dllexport). On ELF/Mach-O both expand to nothing useful --
# FL_API is visibility("default") regardless -- so this costs nothing there.
src/capi.o:     FL_API_DEF = -DFLATLAND_STATIC=1
src/capi.pic.o: FL_API_DEF = -DFLATLAND_BUILD_SHARED=1

-include $(DEP)

test: $(BIN) $(STATICLIB)
	./tests/run_tests.sh ./$(BIN)

# The full validation study against closed-form results, including the
# convergence sweeps behind docs/VALIDATION.md. `make test` runs a quick subset.
validate: $(BIN)
	$(PYTHON) validation/self_check.py
	$(PYTHON) validation/run_validation.py --flatland ./$(BIN)

install: $(BIN) lib
	install -d $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/
	install -m 644 $(STATICLIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(SHAREDLIB) $(DESTDIR)$(PREFIX)/lib/
	$(if $(IMPLIB),install -m 644 $(IMPLIB) $(DESTDIR)$(PREFIX)/lib/)
	install -m 644 include/flatland.h $(DESTDIR)$(PREFIX)/include/

clean:
	rm -f $(BIN) $(STATICLIB) $(SHAREDLIB) $(IMPLIB) \
	      $(CORE_OBJ) $(CAPI_OBJ) $(CLI_OBJ) $(LIB_PIC) $(DEP) \
	      tests/capi_test
