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
SRC       = src/main.cpp
BIN       = flatland

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(BIN) $(LDFLAGS)

.PHONY: test clean
test: $(BIN)
	./tests/run_tests.sh ./$(BIN)

clean:
	rm -f $(BIN)
