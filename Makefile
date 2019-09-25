SHELL := bash
CXX := g++
CPPFLAGS := -std=c++14 -Iinclude
CXXFLAGS := -Wall -O3 -flto -fmax-errors=3 $(CPPFLAGS)
# CXXFLAGS := -Wall -g -fmax-errors=3 $(CPPFLAGS) -DDEBUG_AT
LDFLAGS :=
LDLIBS :=

BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS  := $(shell root-config --ldflags)
ROOT_LDLIBS   := $(shell root-config --libs)

SRCS := $(shell find -L src -type f -name '*$(EXT)')
DEPS := $(patsubst src/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' src --include='*$(EXT)'
EXES := $(patsubst src%$(EXT),bin%, $(shell $(GREP_EXES)))

all: $(EXES)

-include $(DEPS)

# -------------------------------------------------------------------
C_mxaod_4vec := $(ROOT_CXXFLAGS)
L_mxaod_4vec := $(ROOT_LDLIBS) -lTreePlayer
C_varcmp_data := $(ROOT_CXXFLAGS)
L_varcmp_data := $(ROOT_LDLIBS) -lTreePlayer
C_varcmp_mc := $(ROOT_CXXFLAGS)
L_varcmp_mc := $(ROOT_LDLIBS) -lTreePlayer
C_mxaod_4vec2 := $(ROOT_CXXFLAGS)
L_mxaod_4vec2 := $(ROOT_LDLIBS) -lTreePlayer -lpcre

bin/read2 bin/filter2 bin/bin2: \
  $(BLD)/ivanp/io/mem_file.o
# -------------------------------------------------------------------

$(DEPS): $(BLD)/%.d: src/%$(EXT)
	@mkdir -pv $(dir $@)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

bin/%: $(BLD)/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

endif

clean:
	@rm -rfv $(BLD) bin

