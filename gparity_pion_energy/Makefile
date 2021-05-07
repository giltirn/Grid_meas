SHELL := /bin/bash
PREFIX := $(shell grid-config --prefix)
CXX := $(shell grid-config --cxx)
CXXLD := $(shell grid-config --cxxld)
CXXFLAGS := $(shell grid-config --cxxflags)
LDFLAGS := $(shell grid-config --ldflags)
LIBS := $(shell grid-config --libs) $(PREFIX)/lib/libGrid.a
INCLUDES := -I$(PREFIX)/include


main: main.cc
	$(CXX) main.cc -c -o main.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) main.o -o main $(CXXFLAGS) $(LDFLAGS) $(LIBS)
