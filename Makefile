CFLAGS = -std=c++2a -Wall -Wextra -Wshadow -pedantic -march=native

PERF_FLAGS = -Ofast -DNDEBUG -fopenmp

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE = -isystem deps/sdsl-lite/include -isystem deps/seqio/include

HEADERS = include/gapmer.hpp include/fm_index.hpp

.PHONY: clean all fast debug

.DEFAULT: all

%/%.hpp:

all: fast debug

fast: huddinge

debug: huddinge_deb

huddinge: huddinge.cpp $(HEADERS)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge

huddinge_deb: huddinge.cpp $(HEADERS)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge_deb