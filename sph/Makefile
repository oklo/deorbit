CXX ?= clang++
CXXFLAGS ?= -O2 -std=c++17 -Wall -pthread

HDRS = sph.hpp vec3.hpp eos.hpp

sod: sod.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) sod.cpp -o sod

eos_test: eos_test.cpp eos.hpp
	$(CXX) $(CXXFLAGS) eos_test.cpp -o eos_test

.PHONY: clean
clean:
	rm -f sod eos_test sod_out.txt
