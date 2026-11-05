# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Jan Suchánek, Baremon s.r.o.
#
# make                  library, tests, benchmarks and example in $(BUILD)
# make test             correctness tests, a few minutes (AVX-512 IFMA, or AVX-512F with EMULATE_IFMA=1)
# make test-full        the same, test_bpsw also with the timing of Table 3
#
# Variables: CXX/CC (default clang++/clang, as in the paper), ARCH (default -march=native), SQUARING=0 disables the
# dedicated squaring of the filter, GMP_INCLUDE (directory of gmp.h if it is not in the default path),
# GMP_LIB (default -lgmp), CRYPTO_LIB (default -lcrypto), BUILD (default build).
# GMP_RPATH optionally records the directory containing the selected runtime libgmp.
# EMULATE_IFMA=1 replaces the two IFMA instructions by an exact AVX-512F emulation, so that the tests run on CPUs
# without IFMA (use a separate BUILD directory; the timings of such a build mean nothing).

# make predefines CXX and CC (g++, cc), so ?= would not take effect for them
ifeq ($(origin CXX),default)
CXX := clang++
endif
ifeq ($(origin CC),default)
CC := clang
endif
ARCH       ?= -march=native
OPT        ?= -O3
SQUARING   ?= 1
BUILD      ?= build
GMP_LIB    ?= -lgmp
GMP_RPATH  ?=
CRYPTO_LIB ?= -lcrypto
comma := ,
GMP_RPATH_FLAG = $(if $(GMP_RPATH),-Wl$(comma)-rpath$(comma)$(GMP_RPATH))

EMULATE_IFMA ?= 0

CPPFLAGS += -Iinclude $(if $(GMP_INCLUDE),-I$(GMP_INCLUDE)) $(if $(filter 1,$(SQUARING)),-DUSE_AMS) \
            $(if $(filter 1,$(EMULATE_IFMA)),-DLANEBPSW_EMULATE_IFMA)
CXXFLAGS += $(OPT) $(ARCH) -std=c++17
CFLAGS   += $(OPT) $(ARCH)

HEADERS := $(wildcard include/lanebpsw/*.hpp) include/lanebpsw.h
TESTS   := test_kernels test_bpsw test_lucv test_paths test_hashprime test_api
BENCHES := bench_np3 bench_pass prof_np2 cyc_ams bench_sha bench_gmp_builds

.PHONY: all lib tests benches test test-full clean FORCE
all: lib tests benches $(BUILD)/example
lib: $(BUILD)/liblanebpsw.so
tests: $(addprefix $(BUILD)/,$(TESTS))
benches: $(addprefix $(BUILD)/,$(BENCHES))

$(BUILD):
	mkdir -p $@

FORCE:
$(BUILD)/.config: FORCE | $(BUILD)
	@{ printf '%s\n' \
	  'CC=$(CC)' 'CC_VERSION=$(shell $(CC) --version 2>/dev/null | head -1)' \
	  'CXX=$(CXX)' 'CXX_VERSION=$(shell $(CXX) --version 2>/dev/null | head -1)' \
	  'CPPFLAGS=$(CPPFLAGS)' 'CFLAGS=$(CFLAGS)' 'CXXFLAGS=$(CXXFLAGS)' \
	  'ARCH=$(ARCH)' 'OPT=$(OPT)' 'SQUARING=$(SQUARING)' 'EMULATE_IFMA=$(EMULATE_IFMA)' \
	  'GMP_INCLUDE=$(GMP_INCLUDE)' 'GMP_LIB=$(GMP_LIB)' 'GMP_RPATH=$(GMP_RPATH)' \
	  'CRYPTO_LIB=$(CRYPTO_LIB)'; } > $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@
	@rm -f $@.tmp

$(BUILD)/cpu.o: src/cpu.c include/lanebpsw.h $(BUILD)/.config | $(BUILD)
	$(CC) $(CPPFLAGS) $(OPT) -fPIC -c -o $@ $<

$(BUILD)/lanebpsw.o: src/lanebpsw.cpp $(HEADERS) $(BUILD)/.config | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fPIC -c -o $@ $<

$(BUILD)/liblanebpsw.so: $(BUILD)/cpu.o $(BUILD)/lanebpsw.o | $(BUILD)
	$(CXX) -shared -o $@ $^ $(GMP_LIB) $(CRYPTO_LIB) $(GMP_RPATH_FLAG)

$(BUILD)/test_api: tests/test_api.c $(BUILD)/liblanebpsw.so $(BUILD)/.config | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< -L$(BUILD) -llanebpsw -Wl,-rpath,'$$ORIGIN' $(GMP_LIB) $(GMP_RPATH_FLAG)

$(BUILD)/test_hashprime: tests/test_hashprime.cpp $(HEADERS) $(BUILD)/.config | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(GMP_LIB) $(CRYPTO_LIB) $(GMP_RPATH_FLAG)

$(BUILD)/example: examples/example.c $(BUILD)/liblanebpsw.so $(BUILD)/.config | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< -L$(BUILD) -llanebpsw -Wl,-rpath,'$$ORIGIN' $(GMP_LIB) $(GMP_RPATH_FLAG)

$(BUILD)/test_%: tests/test_%.cpp $(HEADERS) $(BUILD)/.config | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(GMP_LIB) $(GMP_RPATH_FLAG)

$(BUILD)/bench_sha: bench/bench_sha.c $(BUILD)/.config | $(BUILD)
	$(CC) $(CFLAGS) -o $@ $< $(CRYPTO_LIB)

$(BUILD)/bench_gmp_builds: bench/bench_gmp_builds.c $(BUILD)/.config | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< -ldl

$(BUILD)/%: bench/%.cpp $(HEADERS) $(BUILD)/.config | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(GMP_LIB) $(GMP_RPATH_FLAG)

test: tests
	python3 tests/test_hp_measurement_gate.py
	$(BUILD)/test_kernels
	$(BUILD)/test_api
	$(BUILD)/test_lucv
	$(BUILD)/test_paths
	$(BUILD)/test_hashprime
	$(BUILD)/test_bpsw --checks-only

test-full: tests
	python3 tests/test_hp_measurement_gate.py
	$(BUILD)/test_kernels
	$(BUILD)/test_api
	$(BUILD)/test_lucv
	$(BUILD)/test_paths
	$(BUILD)/test_hashprime
	$(BUILD)/test_bpsw

clean:
	rm -rf $(BUILD)
