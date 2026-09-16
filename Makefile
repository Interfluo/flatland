# FlatLand — dependency-free build.
#
# No CMake, no third-party libraries: just a C++17 compiler. Ideal for airgapped
# or minimal systems.
#
#   make            # build ./flatland
#   make test       # build, then run the test suite
#   make clean      # remove build artifacts
#
# Override the compiler if needed:  make CXX=clang++

CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall
LDFLAGS  ?= -pthread
INCLUDES  = -Isrc
SRC       = src/fl_mesh.cpp \
            src/fl_io.cpp \
            src/fl_field.cpp \
            src/fl_raster.cpp \
            src/fl_project.cpp \
            src/fl_image.cpp \
            src/fl_batch.cpp \
            cli/main.cpp
OBJ       = $(SRC:.cpp=.o)
DEP       = $(OBJ:.o=.d)
BIN       = flatland

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(BIN) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

-include $(DEP)

.PHONY: test clean
test: $(BIN)
	./tests/run_tests.sh ./$(BIN)

clean:
	rm -f $(BIN) $(OBJ) $(DEP)
