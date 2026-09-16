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
PREFIX   ?= /usr/local

INCLUDES  = -Isrc -Iinclude

# The engine. Shared by the CLI, the C API and anything else that links it.
CORE_SRC = src/fl_mesh.cpp \
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

BIN       = flatland
STATICLIB = libflatland.a

# macOS wants .dylib and a different soname flag; everything else gets .so.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHAREDLIB = libflatland.dylib
  SOFLAGS   = -dynamiclib -install_name @rpath/$(SHAREDLIB)
else
  SHAREDLIB = libflatland.so
  SOFLAGS   = -shared -Wl,-soname,$(SHAREDLIB)
endif

.PHONY: all lib test test-cli test-capi install clean

$(BIN): $(CORE_OBJ) $(CLI_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

all: $(BIN) lib
lib: $(STATICLIB) $(SHAREDLIB)

$(STATICLIB): $(CORE_OBJ) $(CAPI_OBJ)
	$(AR) rcs $@ $^

# -fvisibility=hidden keeps the engine's C++ symbols internal, so the shared
# library exports only what flatland.h marks FL_API.
$(SHAREDLIB): $(LIB_PIC)
	$(CXX) $(CXXFLAGS) $(SOFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

%.pic.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -fPIC -fvisibility=hidden -MMD -MP -c $< -o $@

-include $(DEP)

test: $(BIN) $(STATICLIB)
	./tests/run_tests.sh ./$(BIN)

install: $(BIN) lib
	install -d $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/include
	install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/
	install -m 644 $(STATICLIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(SHAREDLIB) $(DESTDIR)$(PREFIX)/lib/
	install -m 644 include/flatland.h $(DESTDIR)$(PREFIX)/include/

clean:
	rm -f $(BIN) $(STATICLIB) $(SHAREDLIB) \
	      $(CORE_OBJ) $(CAPI_OBJ) $(CLI_OBJ) $(LIB_PIC) $(DEP) \
	      tests/capi_test
