-include local.mk

ifndef MAX_GAP
MAX_GAP = 15
endif

GCOV ?= gcov

CFLAGS = -std=c++23 -Wall -Wextra -Wshadow -pedantic -march=native -DMAX_GAP=$(MAX_GAP)

PERF_FLAGS = -Ofast -DNDEBUG -fopenmp
TEST_PERF_FLAGS = -O0 -g -DDEBUG

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE += -isystem deps/sdsl-lite/include -isystem deps/seqio/include

LIBS += -lz -lgsl -lgslcblas -lm

HEADERS = include/gapmer.hpp include/fm_index.hpp include/gapmer_count.hpp include/seed_finder.hpp include/partial_count.hpp include/seed_clusterer.hpp

SDSL_DIR = deps/sdsl-lite/lib

GTEST_DIR = deps/googletest
RAPIDCHECK_DIR = deps/rapidcheck

GFLAGS = -lpthread -DGTEST_ON -isystem $(GTEST_DIR)/googletest/include -isystem $(RAPIDCHECK_DIR)/include -isystem $(RAPIDCHECK_DIR)/extras/gtest/include -pthread -L $(GTEST_DIR)/build/lib -L $(RAPIDCHECK_DIR)/build

TEST_HPP = test/gapmer_tests.hpp test/util_tests.hpp test/gapmer_count_tests.hpp

.PHONY: clean all fast debug test cover

.DEFAULT: all

%/%.hpp:

all: fast debug

all_: huddinge seed_finder comp huddinge_deb seed_finder_deb

fast: huddinge

debug: huddinge_deb

huddinge: huddinge.cpp include/util.hpp include/gapmer.hpp
	$(CXX) $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge

seed_finder: seed_finder.cpp $(HEADERS) | $(SDSL_DIR)
	$(CXX) $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) seed_finder.cpp -o seed_finder $(LIBS)

comp: comp.cpp include/util.hpp
	$(CXX) $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) comp.cpp -o comp

huddinge_deb: huddinge.cpp $(HEADERS)
	$(CXX) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge_deb

seed_finder_deb: seed_finder.cpp $(HEADERS) | $(SDSL_DIR)
	$(CXX) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) seed_finder.cpp -o seed_finder_deb $(LIBS)

clean:
	rm -f huddinge huddinge_deb seed_finder seed_finder_deb comp
	rm -f test/test test/cover test/test.o test/cover.o
	rm -f *.gcov test/*.gcda test/*.gcno index.info
	rm -rf target/

$(SDSL_DIR):
	git submodule update --init

$(GTEST_DIR)/googletest:
	git submodule update --init

$(RAPIDCHECK_DIR):
	git submodule update --init

$(GTEST_DIR)/build/lib/libgtest_main.a: | $(GTEST_DIR)/googletest
	(mkdir -p $(GTEST_DIR)/build && cd $(GTEST_DIR)/build && cmake -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF .. && $(MAKE))

$(RAPIDCHECK_DIR)/build/librapidcheck.a: | $(RAPIDCHECK_DIR)
	(mkdir -p $(RAPIDCHECK_DIR)/build && cd $(RAPIDCHECK_DIR)/build && cmake -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF -DRC_ENABLE_GTEST=ON .. && $(MAKE))

test/test.o: $(TEST_HPP) $(GTEST_HEADERS) $(HEADERS) test/test.cpp
	$(CXX) $(CFLAGS) $(GFLAGS) -c test/test.cpp -o test/test.o

test/test: $(GTEST_DIR)/build/lib/libgtest_main.a $(RAPIDCHECK_DIR)/build/librapidcheck.a $(TEST_HPP) $(HEADERS) test/test.cpp
	$(CXX) $(CFLAGS) $(GFLAGS) $(INCLUDE) $(TEST_PERF_FLAGS) test/test.cpp -o test/test -lgtest_main -lgtest -lrapidcheck $(LIBS)

test: test/test
	test/test $(ARG)

test/cover: $(GTEST_DIR)/build/lib/libgtest_main.a $(TEST_HPP) $(HEADERS) test/test.cpp
	$(CXX) -g --coverage -O1 $(CFLAGS) $(GFLAGS) $(INCLUDE) test/test.cpp -o test/cover -lgtest_main -lgtest -lrapidcheck -lgcov $(LIBS)

cover: test/cover
	rm -f *.gcov test/*.gcda coverage.info
	test/cover
	lcov -c --ignore-errors inconsistent,unused --directory test --output-file coverage.info --gcov-tool $(GCOV)
	lcov --remove coverage.info "/usr*" --output-file coverage.info
	genhtml coverage.info -o target
