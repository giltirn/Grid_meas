SHELL := /bin/bash

SOURCE := /direct/sdcc+u/ckelly/MLcluster/src/Grid_meas/Grid_meas

PREFIX := $(shell grid-config --prefix)
CXX := $(shell grid-config --cxx)
CXXLD := $(shell grid-config --cxxld)
CXXFLAGS := $(shell grid-config --cxxflags) -g
LDFLAGS := $(shell grid-config --ldflags) -g
LIBS := $(shell grid-config --libs) $(PREFIX)/lib/libGrid.a
INCLUDES := -I$(PREFIX)/include -I$(SOURCE)

test_sources: test_sources.cc
	$(CXX) test_sources.cc -c -o test_sources.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_sources.o -o test_sources $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_propagator_invert: test_propagator_invert.cc
	$(CXX) test_propagator_invert.cc -c -o test_propagator_invert.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_propagator_invert.o -o test_propagator_invert $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_gparity_cconj_reln: test_gparity_cconj_reln.cc
	$(CXX) test_gparity_cconj_reln.cc -c -o test_gparity_cconj_reln.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_gparity_cconj_reln.o -o test_gparity_cconj_reln $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_momentum: test_momentum.cc
	$(CXX) test_momentum.cc -c -o test_momentum.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_momentum.o -o test_momentum $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_array_io: test_array_io.cc
	$(CXX) test_array_io.cc -c -o test_array_io.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_array_io.o -o test_array_io $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_pion_2pt: test_pion_2pt.cc
	$(CXX) test_pion_2pt.cc -c -o test_pion_2pt.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_pion_2pt.o -o test_pion_2pt $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_kaon_2pt: test_kaon_2pt.cc
	$(CXX) test_kaon_2pt.cc -c -o test_kaon_2pt.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_kaon_2pt.o -o test_kaon_2pt $(CXXFLAGS) $(LDFLAGS) $(LIBS)

test_split_grid_propagator_invert: test_split_grid_propagator_invert.cc
	$(CXX) test_split_grid_propagator_invert.cc -c -o test_split_grid_propagator_invert.o $(INCLUDES) $(CXXFLAGS)
	$(CXXLD) test_split_grid_propagator_invert.o -o test_split_grid_propagator_invert $(CXXFLAGS) $(LDFLAGS) $(LIBS)
