-include local.mk

ifndef MAX_GAP
MAX_GAP = 15
endif

GCOV ?= gcov

WARNING_FLAGS = -Wall -Werror -Wextra -Wshadow -Wno-gnu-conditional-omitted-operand -Wno-unused-parameter -Wno-unused-function -Wno-error=unused-command-line-argument
CFLAGS = -std=c++23 $(WARNING_FLAGS) -march=native -DMAX_GAP=$(MAX_GAP)

PERF_FLAGS = -O3 -g -DNDEBUG -fopenmp
TEST_PERF_FLAGS = -O0 -g -DDEBUG
TEST_COVERAGE_PERF_FLAGS = -O1 -g --coverage

DEBUG_FLAGS = -g -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE += -Iinclude -isystem deps/sdsl-lite/include -isystem deps/seqio/include -isystem deps/args -isystem deps/libbio/include -isystem deps/libbio/lib/range-v3/include

LIBS += -lboost_iostreams -lz -lgsl -lgslcblas -lm

HEADERS =	include/args.hpp \
			include/bits.hpp \
			include/configuration.hpp \
			include/dot_writer.hpp \
			include/fm_index.hpp \
			include/gapmer.hpp \
			include/gapmer_count.hpp \
			include/libbio_reader_adapter.hpp \
			include/pack_characters.hpp \
			include/packed_character_iteration.hpp \
			include/partial_count.hpp \
			include/reader_adapter.hpp \
			include/seed_clusterer.hpp \
			include/seed_finder.hpp \
			include/seqio_reader_adapter.hpp \
			include/string_buffer.hpp \
			include/util.hpp \
			include/version.hpp

# FIXME: A library-based design could be a good idea.
SEED_FINDER_OBJECTS = \
			pack_characters.o \
			seed_finder.o \
			seqio_reader_adapter.o

TEST_OBJECTS = \
			pack_characters.o \
			test/test.o

TEST_COVERAGE_OBJECTS = $(TEST_OBJECTS:.o=.coverage.o)
GCDA = $(TEST_COVERAGE_OBJECTS:.o=.gcda)
GCNO = $(TEST_COVERAGE_OBJECTS:.o=.gcno)

ARGS = deps/args/args.hxx
LIBBIO_DIR = deps/libbio
SDSL_DIR = deps/sdsl-lite/lib

GTEST_DIR = deps/googletest
RAPIDCHECK_DIR = deps/rapidcheck

TEST_CPPFLAGS = -DGTEST_ON -isystem $(GTEST_DIR)/googletest/include -isystem $(RAPIDCHECK_DIR)/include -isystem $(RAPIDCHECK_DIR)/extras/gtest/include
TEST_LDFLAGS = -pthread -L $(GTEST_DIR)/build/lib -L $(RAPIDCHECK_DIR)/build

TEST_HPP =	test/bit_tests_arbitrary.hpp \
			test/gapmer_arbitrary.hpp \
			test/gapmer_count_tests.hpp \
			test/gapmer_tests.hpp \
			test/nucleotide.hpp \
			test/pack_characters_arbitrary.hpp \
			test/packed_character_iteration_arbitrary.hpp \
			test/string_buffer_arbitrary.hpp \
			test/test.hpp \
			test/util_tests.hpp


.PHONY: clean all test cover

%/%.hpp:

motivating: motivating.cpp
	$(CXX) $(CFLAGS) $(PERF_FLAGS) -isystem deps/seqio/include motivating.cpp -o motivating

huddinge: huddinge.cpp include/util.hpp include/gapmer.hpp
	$(CXX) $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge

seed_finder: $(SEED_FINDER_OBJECTS) $(HEADERS) $(ARGS) | $(SDSL_DIR) $(ARGS_DIR) $(LIBBIO_DIR)
	$(CXX) $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) $(SEED_FINDER_OBJECTS) -o seed_finder $(LIBS)

comp: comp.cpp include/util.hpp
	$(CXX) $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) comp.cpp -o comp

huddinge_deb: huddinge.cpp $(HEADERS)
	$(CXX) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) huddinge.cpp -o huddinge_deb

seed_finder_deb: $(SEED_FINDER_OBJECTS) $(HEADERS) $(ARGS) | $(SDSL_DIR) $(ARGS_DIR) $(LIBBIO_DIR)
	$(CXX) $(CFLAGS) $(DEBUG_FLAGS) $(INCLUDE) $(SEED_FINDER_OBJECTS) -o seed_finder_deb $(LIBS)

clean:
	$(RM) $(SEED_FINDER_OBJECTS) $(TEST_OBJECTS) $(TEST_COVERAGE_OBJECTS)
	$(RM) $(GCDA) $(GCNO) coverage.info
	rm -f huddinge huddinge_deb seed_finder seed_finder_deb comp
	rm -f test/test test/cover
	rm -rf target/

$(SDSL_DIR):
	git submodule update --init

$(GTEST_DIR)/googletest:
	git submodule update --init

$(RAPIDCHECK_DIR):
	git submodule update --init

$(ARGS):
	git submodule update --init

$(LIBBIO_DIR):
	git submodule update --init --recursive

$(LIBBIO_DIR)/src/libbio.a: deps/libbio/local.mk
	$(MAKE) -C deps/libbio src/libbio.a

$(LIBBIO_DIR)/local.mk: local.mk
	cp local.mk deps/libbio/

local.mk:
	bash -c "[ -e local.mk ] || echo -n "" > local.mk"

$(GTEST_DIR)/build/lib/libgtest_main.a: | $(GTEST_DIR)/googletest
	(mkdir -p $(GTEST_DIR)/build && cd $(GTEST_DIR)/build && cmake -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF .. && $(MAKE))

$(RAPIDCHECK_DIR)/build/librapidcheck.a: | $(RAPIDCHECK_DIR)
	(mkdir -p $(RAPIDCHECK_DIR)/build && cd $(RAPIDCHECK_DIR)/build && cmake -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF -DRC_ENABLE_GTEST=ON .. && $(MAKE))

test: test/test
	test/test $(ARG)

test/test: $(GTEST_DIR)/build/lib/libgtest_main.a $(RAPIDCHECK_DIR)/build/librapidcheck.a $(HEADERS) $(TEST_HPP) $(TEST_OBJECTS)
	$(CXX) $(CFLAGS) $(INCLUDE) $(TEST_PERF_FLAGS) $(TEST_OBJECTS) -o test/test $(TEST_LDFLAGS) -lgtest_main -lgtest -lrapidcheck $(LIBS)

test/cover: $(GTEST_DIR)/build/lib/libgtest_main.a $(RAPIDCHECK_DIR)/build/librapidcheck.a $(HEADERS) $(TEST_HPP) $(TEST_COVERAGE_OBJECTS)
	$(CXX) $(CFLAGS) $(INCLUDE) $(TEST_COVERAGE_PERF_FLAGS) $(TEST_COVERAGE_OBJECTS) -o test/cover $(TEST_LDFLAGS) -lgtest_main -lgtest -lrapidcheck -lgcov $(LIBS)

cover: test/cover
	$(RM) $(GCDA) $(GCNO) coverage.info
	test/cover
	lcov -c --ignore-errors inconsistent,unused --directory test --output-file coverage.info --gcov-tool $(GCOV)
	lcov --remove coverage.info "/usr*" --output-file coverage.info
	genhtml coverage.info -o target

test/test.o: $(TEST_HPP) $(GTEST_HEADERS) $(HEADERS) test/test.cpp
	$(CXX) $(CFLAGS) $(TEST_CPPFLAGS) $(TEST_PERF_FLAGS) $(INCLUDE) -c test/test.cpp -o test/test.o

test/test.coverage.o: $(TEST_HPP) $(GTEST_HEADERS) $(HEADERS) test/test.cpp
	$(CXX) $(CFLAGS) $(TEST_CPPFLAGS) $(TEST_COVERAGE_PERF_FLAGS) $(INCLUDE) -c test/test.cpp -o test/test.coverage.o

%.o: %.cpp
	$(CXX) -c -o $@ $(CFLAGS) $(PERF_FLAGS) $(INCLUDE) $<

%.coverage.o: %.cpp
	$(CXX) -c -o $@ $(CFLAGS) $(TEST_COVERAGE_PERF_FLAGS) $(INCLUDE) $<

# Output git tag or version hash to config.h.
include/version.hpp: SHELL := /bin/bash
include/version.hpp: .git
	VERSION_STR=`./git_version.sh`; \
	printf "// This file is automatically generated by Makefile.\n#include <string>\nnamespace sf {\n  constexpr static std::string const version{\"%s\"};\n}\n" "$${VERSION_STR}" > "$@"
