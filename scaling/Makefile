# File       : Makefile
# Description: Makefile for peak performance code
# Copyright 2023 Harvard University. All Rights Reserved.
CXX = g++
CPPFLAGS_DEBUG = -g -O0 -fsanitize=address -Wall -Wextra -fopenmp
CPPFLAGS = -O3 -fopenmp -march=native -mtune=native -flto -funroll-loops -ffast-math

.PHONY: clean

strong: StrongScaling.cpp
	$(CXX) $(CPPFLAGS) -o $@ $^

weak: WeakScaling.cpp
	$(CXX) $(CPPFLAGS) -o $@ $^
