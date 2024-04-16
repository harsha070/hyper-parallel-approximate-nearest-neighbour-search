# File       : Makefile
# Description: Makefile for peak performance code
# Copyright 2023 Harvard University. All Rights Reserved.
CXX = g++
MPICXX = mpic++
# CXXFLAGS = -g -O0 -fsanitize=address -Wall -Wextra -fopenmp #
CXXFLAGS = -g -O0 -Wall -Wextra -fopenmp -gdwarf-2 -gstrict-dwarf -gdwarf-4
CPPFLAGS_FAST = -O3 -fopenmp


# TODO: add compiler flags for your kernel here
# PERFLAGS = -O3 -msse -msse2 -msse3 -DNDEBUG # targets: SSE, SSE2 and SSE3
# PERFLAGS = -O3 -mavx2 -mfma  -ffast-math #-msse3 -march=native -mtune=native -ftree-vectorize
.PHONY: clean

test: test.cpp SearchTree.o
	$(CXX) $(CXXFLAGS) -o $@ $^

SearchTree.o: core/SearchTree.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

main: main.cpp SearchForest.o SearchTree.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^

inference: inference.cpp SearchForest.o SearchTree.o
	$(MPICXX) $(CXXFLAGS) -fsanitize=address -o $@ $^

fast: main.cpp SearchTree.o
	$(CXX) $(CPPFLAGS_FAST) -o $@ $^

SearchForest.o: core/SearchForest.cpp core/SearchForest.hpp
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

forest: main.cpp SearchTree.o SearchForest.o
	$(MPICXX) $(CXXFLAGS) -o $@ $^

clean:
	rm test SearchForest.o SearchTree.o forest
