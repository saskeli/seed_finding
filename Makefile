CFLAGS = -std=c++2a -Wall -Wextra -Wshadow -pedantic -march=native

PERF_FLAGS = -Ofast -DNDEBUG

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

HEADERS = include/kmer.hpp

.PHONY: clean all fast debug

.DEFAULT: all

%/%.hpp:

all: fast debug

fast: huddinge

debug: huddinge_deb

huddinge: huddinge.cpp $(HEADERS)
	g++ $(CFLAGS) $(PERF_FLAGS) huddinge.cpp -o huddinge

huddinge_deb: huddinge.cpp $(HEADERS)
	g++ $(CFLAGS) $(DEBUG_FLAGS) huddinge.cpp -o huddinge_deb