CXX ?= clang++
CXXFLAGS ?= -O2 -std=c++17 -Wall

sod: sod.cpp sph.hpp vec3.hpp
	$(CXX) $(CXXFLAGS) sod.cpp -o sod

.PHONY: clean
clean:
	rm -f sod sod_out.txt
