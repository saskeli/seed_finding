ifndef MAX_GAP
MAX_GAP = 10
endif

CFLAGS = -std=c++23 -Wall -Wextra -Wshadow -pedantic -march=native -DMAX_GAP=$(MAX_GAP)

PERF_FLAGS = -Ofast -DNDEBUG -fopenmp

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE = -isystem deps/sdsl-lite/include -isystem deps/seqio/include

LIBS = -L deps/sdsl-lite/lib -lsdsl -lz -lgsl -lgslcblas -lm

HEADERS = include/gapmer.hpp include/fm_index.hpp include/gapmer_count.hpp include/seed_finder.hpp

SDSL_DIR = deps/sdsl-lite/lib

SDSL_A = $(SDSL_DIR)/libsdsl.a

GTEST_DIR = deps/googletest

GFLAGS = -lpthread -DGTEST_ON -isystem $(GTEST_DIR)/googletest/include -pthread -L $(GTEST_DIR)/lib

GTEST_HEADERS = $(GTEST_DIR)/googletest/include/gtest/*.h \
                $(GTEST_DIR)/googletest/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/googletest/src/*.cc $(GTEST_DIR)/googletest/src/*.h $(GTEST_HEADERS)

TEST_HPP = test/gapmer_tests.hpp test/util_tests.hpp test/gapmer_count_tests.hpp

.PHONY: clean all fast debug test cover

.DEFAULT: all

%/%.hpp:

all: fast debug

fast: huddinge

debug: huddinge_deb

huddinge: huddinge.cpp include/util.hpp include/gapmer.hpp
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge

seed_finder: seed_finder.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) seed_finder.cpp -o seed_finder $(LIBS)

comp: comp.cpp include/util.hpp
	g++ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) comp.cpp -o comp

huddinge_deb: huddinge.cpp $(HEADERS)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge_deb

seed_finder_deb: seed_finder.cpp $(HEADERS) $(SDSL_A)
	g++ $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) seed_finder.cpp -o seed_finder_deb $(LIBS)

clean:
	rm -f huddinge huddinge_deb seed_finder seed_finder_deb
	rm -f test/test test/cover test/test.o test/cover.o
	rm -f *.gcov test/*.gcda test/*.gcno index.info
	rm -rf target/

$(SDSL_DIR):
	git submodule init
	git submodule update

$(GTEST_DIR)/googletest:
	git submodule init
	git submodule update

$(GETST_DIR)/lib/libgtest_main.a: | $(GTEST_DIR)/googletest
	(cd $(GTEST_DIR) && cmake CMakelists.txt && make)

$(SDSL_A): | $(SDSL_DIR)
	(cd deps/sdsl-lite && cmake CMakelists.txt && make)

test/test.o: $(TEST_HPP) $(GTEST_HEADERS) $(HEADERS) test/test.cpp
	g++ $(CFLAGS) $(GFLAGS) -c test/test.cpp -o test/test.o

test/test: $(GTEST_DIR)/lib/libgtest_main.a $(TEST_HPP) $(GTEST_HEADERS) $(HEADERS) test/test.cpp
	g++ $(CFLAGS) $(GFLAGS) $(INCLUDE) -O3 test/test.cpp -o test/test -lgtest_main -lgtest $(LIBS)

test: test/test
	test/test $(ARG)

test/cover: $(GTEST_DIR)/lib/libgtest_main.a $(TEST_HPP) $(GTEST_HEADERS) $(HEADERS) test/test.cpp
	g++ -g --coverage -O1 $(CFLAGS) $(GFLAGS) test/test.cpp -o test/cover -lgtest_main -lgtest -lgcov

cover: test/cover
	rm -f *.gcov test/*.gcda coverage.info
	test/cover
	lcov -c --ignore-errors inconsistent,unused --directory test --output-file coverage.info --gcov-tool gcov
	lcov --remove coverage.info "/usr*" --output-file coverage.info
	genhtml coverage.info -o target

