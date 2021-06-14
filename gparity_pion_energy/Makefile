SHELL := /bin/bash

SOURCE := /direct/sdcc+u/ckelly/MLcluster/src/Grid_meas/Grid_meas

PREFIX := $(shell grid-config --prefix)
CXX := $(shell grid-config --cxx)
CXXLD := $(shell grid-config --cxxld)
CXXFLAGS := $(shell grid-config --cxxflags) -g
LDFLAGS := $(shell grid-config --ldflags) -g
LIBS := $(shell grid-config --libs) $(PREFIX)/lib/libGrid.a
INCLUDES := -I$(PREFIX)/include -I$(SOURCE)

main: main.cc
	$(CXX) main.cc -c -o main.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) main.o -o main $(CXXFLAGS) $(LDFLAGS) $(LIBS)
