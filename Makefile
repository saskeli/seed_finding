CFLAGS = -std=c++2a -Wall -Wextra -Wshadow -pedantic -march=native

PERF_FLAGS = -Ofast -DNDEBUG -fopenmp

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE = -isystem deps/sdsl-lite/include -isystem deps/seqio/include

LIBS = -L deps/sdsl-lite/lib -lsdsl -lz

HEADERS = include/gapmer.hpp include/fm_index.hpp

SDSL_A = deps/sdsl-lite/lib/libsdsl.a

.PHONY: clean all fast debug

.DEFAULT: all

%/%.hpp:

all: fast debug

fast: huddinge

debug: huddinge_deb

huddinge: huddinge.cpp include/util.hpp
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge

seed_finder: seed_finder.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) seed_finder.cpp -o seed_finder $(LIBS)

comp: comp.cpp include/util.hpp
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) comp.cpp -o comp

huddinge_deb: huddinge.cpp $(HEADERS)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge_deb

seed_finder_deb: seed_finder.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) seed_finder.cpp -o seed_finder_deb $(LIBS)

$(SDSL_A):
	(cd deps/sdsl-lite && cmake CMakelists.txt && make)