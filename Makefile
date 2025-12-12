-include local.mk

ifndef MAX_GAP
MAX_GAP = 15
endif

CMAKE ?= cmake
GCOV ?= gcov

## Not used directly
WARNING_FLAGS_ ?=
WARNING_FLAGS = \
	-Wall \
	-Werror \
	-Wextra \
	-Wshadow \
	-Wno-gnu-conditional-omitted-operand \
	-Wno-unused-parameter \
	-Wno-unused-function \
	$(WARNING_FLAGS_)

PERF_FLAGS = -O3 -g -DNDEBUG -fopenmp -march=native
DEBUG_PERF_FLAGS = -O0 -g -DDEBUG
TEST_PERF_FLAGS = $(DEBUG_PERF_FLAGS)
TEST_COVERAGE_PERF_FLAGS = -O1 -g

DEBUG_FLAGS = -DDEBUG

ifdef VERBOSE
DEBUG_FLAGS += -DVERBOSE
endif

INCLUDE += \
	-Iinclude \
	-isystem deps/sdsl-lite/include \
	-isystem deps/seqio/include \
	-isystem deps/args \
	-isystem deps/libbio/include \
	-isystem deps/libbio/lib/range-v3/include

LIBS += -lboost_iostreams -lz -lgsl -lgslcblas -lm

## Used directly
CPPFLAGS				= -DMAX_GAP=$(MAX_GAP) $(INCLUDE)
CXXFLAGS				= -std=c++23 $(WARNING_FLAGS) $(PERF_FLAGS)
LDFLAGS					= $(LIBS)
DEBUG_CPPFLAGS			= $(CPPFLAGS) $(DEBUG_FLAGS)
DEBUG_CXXFLAGS			= -std=c++23 $(WARNING_FLAGS) $(DEBUG_PERF_FLAGS)
DEBUG_LDFLAGS			= $(LDFLAGS)
TEST_CPPFLAGS			= -DGTEST_ON $(INCLUDE) -isystem $(GTEST_DIR)/googletest/include -isystem $(RAPIDCHECK_DIR)/include -isystem $(RAPIDCHECK_DIR)/extras/gtest/include
TEST_CXXFLAGS			= -std=c++23 -pthread $(WARNING_FLAGS) $(TEST_PERF_FLAGS)
TEST_LDFLAGS			= -L $(GTEST_DIR)/build/lib -L $(RAPIDCHECK_DIR)/build -lgtest_main -lgtest -lrapidcheck $(LIBS)
TEST_COVERAGE_CPPFLAGS	= $(TEST_CPPFLAGS)
TEST_COVERAGE_CXXFLAGS	= -std=c++23 -pthread --coverage $(WARNING_FLAGS) $(TEST_COVERAGE_PERF_FLAGS)
TEST_COVERAGE_LDFLAGS	= $(TEST_LDFLAGS)


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
			libbio_reader_adapter.o \
			pack_characters.o \
			seed_finder.o

SEED_FINDER_DEBUG_OBJECTS = $(SEED_FINDER_OBJECTS:.o=.debug.o)

TEST_OBJECTS = \
			pack_characters.o \
			test/test.o

TEST_COVERAGE_OBJECTS = $(TEST_OBJECTS:.o=.coverage.o)
GCDA = $(TEST_COVERAGE_OBJECTS:.o=.gcda)
GCNO = $(TEST_COVERAGE_OBJECTS:.o=.gcno)

ARGS_HXX = deps/args/args.hxx
LIBBIO_DIR = deps/libbio
SDSL_DIR = deps/sdsl-lite/lib

GTEST_DIR = deps/googletest
RAPIDCHECK_DIR = deps/rapidcheck

TEST_HEADERS = \
			test/bit_tests_arbitrary.hpp \
			test/gapmer_arbitrary.hpp \
			test/gapmer_count_tests.hpp \
			test/gapmer_tests.hpp \
			test/nucleotide.hpp \
			test/pack_characters_arbitrary.hpp \
			test/packed_character_iteration_arbitrary.hpp \
			test/string_buffer_arbitrary.hpp \
			test/test.hpp \
			test/util_tests.hpp


.PHONY: clean all debug test cover


all: motivating huddinge seed_finder comp

debug: seed_finder_deb huddinge_deb

motivating: motivating.o
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

huddinge: huddinge.o include/util.hpp include/gapmer.hpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

seed_finder: $(SEED_FINDER_OBJECTS) $(HEADERS) $(ARGS_HXX) $(LIBBIO_DIR)/src/libbio.a | $(SDSL_DIR) $(ARGS_DIR) $(LIBBIO_DIR)
	$(CXX) $(CXXFLAGS) $(SEED_FINDER_OBJECTS) -o $@ $(LIBBIO_DIR)/src/libbio.a $(LDFLAGS)

comp: comp.o include/util.hpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

huddinge_deb: huddinge.debug.o $(HEADERS)
	$(CXX) $(DEBUG_CXXFLAGS) $< -o $@ $(DEBUG_LDFLAGS)

seed_finder_deb: $(SEED_FINDER_DEBUG_OBJECTS) $(HEADERS) $(ARGS_HXX) | $(SDSL_DIR) $(ARGS_DIR) $(LIBBIO_DIR)
	$(CXX) $(DEBUG_CXXFLAGS) $(SEED_FINDER_DEBUG_OBJECTS) -o $@ $(DEBUG_LDFLAGS)

clean:
	$(RM) $(SEED_FINDER_OBJECTS) motivating.o huddinge.o comp.o $(SEED_FINDER_DEBUG_OBJECTS) huddinge.debug.o $(TEST_OBJECTS) $(TEST_COVERAGE_OBJECTS)
	$(RM) $(GCDA) $(GCNO) coverage.info
	$(RM) huddinge huddinge_deb seed_finder seed_finder_deb comp motivating
	$(RM) test/test test/cover
	$(RM) -r target/

$(SDSL_DIR):
	git submodule update --init

$(GTEST_DIR)/googletest:
	git submodule update --init

$(RAPIDCHECK_DIR):
	git submodule update --init

$(ARGS_HXX):
	git submodule update --init

$(LIBBIO_DIR):
	git submodule update --init --recursive

$(LIBBIO_DIR)/src/libbio.a: deps/libbio/local.mk
	$(MAKE) -C deps/libbio src/libbio.a

$(LIBBIO_DIR)/local.mk: local.mk
	cp local.mk deps/libbio/

local.mk:
	bash -c "[ -e local.mk ] || echo '# Generated automatically according to Makefile.' > local.mk"

$(GTEST_DIR)/build/lib/libgtest_main.a: | $(GTEST_DIR)/googletest
	(mkdir -p $(GTEST_DIR)/build && cd $(GTEST_DIR)/build && $(CMAKE) -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF .. && $(MAKE))

$(RAPIDCHECK_DIR)/build/librapidcheck.a: | $(RAPIDCHECK_DIR)
	(mkdir -p $(RAPIDCHECK_DIR)/build && cd $(RAPIDCHECK_DIR)/build && $(CMAKE) -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF -DRC_ENABLE_GTEST=ON .. && $(MAKE))

test/test: $(TEST_OBJECTS) $(HEADERS) $(TEST_HEADERS) $(GTEST_DIR)/build/lib/libgtest_main.a $(RAPIDCHECK_DIR)/build/librapidcheck.a
	$(CXX) $(TEST_CXXFLAGS) $(TEST_OBJECTS) -o $@ $(TEST_LDFLAGS)

test/cover: $(TEST_COVERAGE_OBJECTS) $(HEADERS) $(TEST_HEADERS) $(GTEST_DIR)/build/lib/libgtest_main.a $(RAPIDCHECK_DIR)/build/librapidcheck.a
	$(CXX) $(TEST_COVERAGE_CXXFLAGS) $(TEST_COVERAGE_OBJECTS) -o $@ $(TEST_COVERAGE_LDFLAGS)

test: test/test
	test/test $(ARG)

cover: test/cover
	$(RM) $(GCDA) $(GCNO) coverage.info
	./test/cover
	lcov -c --ignore-errors inconsistent,unused --directory test --output-file coverage.info --gcov-tool $(GCOV)
	lcov --remove coverage.info "/usr*" --output-file coverage.info
	genhtml coverage.info -o target

test/test.o: test/test.cpp $(HEADERS) $(TEST_HEADERS) $(GTEST_HEADERS)
	$(CXX) $(TEST_CPPFLAGS) $(TEST_CXXFLAGS) -c $< -o $@

test/test.coverage.o: test/test.cpp $(HEADERS) $(TEST_HEADERS) $(GTEST_HEADERS)
	$(CXX) $(TEST_COVERAGE_CPPFLAGS) $(TEST_COVERAGE_CXXFLAGS) -c $< -o $@

%/%.hpp:

%.o: %.cpp
	$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $<

%.debug.o: %.cpp
	$(CXX) -c -o $@ $(DEBUG_CPPFLAGS) $(DEBUG_CXXFLAGS) $<

%.coverage.o: %.cpp
	$(CXX) -c -o $@ $(TEST_COVERAGE_CPPFLAGS) $(TEST_COVERAGE_CXXFLAGS) $<

# Output git tag or version hash to config.h.
include/version.hpp: SHELL := /bin/bash
include/version.hpp: .git
	VERSION_STR=`./git_version.sh`; \
	printf "// This file is automatically generated by Makefile.\n#include <string>\nnamespace sf {\n  constexpr static std::string const version{\"%s\"};\n}\n" "$${VERSION_STR}" > "$@"
