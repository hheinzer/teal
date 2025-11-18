# compiler and libraries
CC = clang
MPICC = OMPI_CC=$(CC) MPICH_CC=$(CC) mpicc
LDLIBS = -lm -lhdf5 -lmetis -lparmetis -lunwind -lunwind-x86_64

# compiler flags
CONFIG ?= debug
CFLAGS = -std=c99 -g -Wall -Wextra -Wpedantic -Wshadow -Wwrite-strings -Wcast-qual -Isrc
ifeq ($(CONFIG), debug)
	CFLAGS += -O0 -fno-omit-frame-pointer -fsanitize-trap -fsanitize=address,undefined
endif
ifeq ($(CONFIG), valgrind)
	CFLAGS += -Og -fno-omit-frame-pointer -DVALGRIND
endif
ifneq (,$(filter $(CONFIG), release gprof))
	CFLAGS += -O3 -march=native -flto=auto -DNDEBUG
endif
ifeq ($(CONFIG), gprof)
	CFLAGS += -pg -fno-inline-functions -fno-optimize-sibling-calls
endif
CFLAGS += -DCONFIG=\"$(CONFIG)\"

# gcc/tcc specific flags (silence bogous warnings)
ifneq (,$(filter $(CC), gcc tcc))
	CFLAGS += -Wno-pedantic -Wno-discarded-qualifiers -Wno-sign-compare
endif

# build info
COMPILER := $(shell $(CC) --version | head -n1)
COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || echo "unknown")
STATUS := $(shell git diff-index --quiet HEAD 2>/dev/null || echo "-dirty")
CFLAGS += -DCOMPILER="\"$(COMPILER)\"" -DCOMMIT=\"$(COMMIT)$(STATUS)\"

# sources, objects, and programs
SRC := $(shell find src -type f -name '*.c')
RUN := $(shell find run -type f -name '*.c')
OBJ := $(patsubst src/%.c, obj/%.o, $(SRC))
BIN := $(patsubst run/%.c, bin/%, $(RUN))

# make functions
.PHONY: all valgrind release gprof clean check tidy format

all: $(BIN)

valgrind:
	@$(MAKE) --no-print-directory CONFIG=valgrind

release:
	@$(MAKE) --no-print-directory CONFIG=release

gprof:
	@$(MAKE) --no-print-directory CONFIG=gprof

clean:
	@rm -rf obj bin

check:
	@cppcheck --quiet --project=compile_commands.json \
		--enable=all --inconclusive --check-level=exhaustive \
		--suppress=checkersReport --suppress=missingIncludeSystem \
		--suppress=constVariable --suppress=constVariablePointer --suppress=unusedFunction

tidy: $(OBJ)
	@clang-tidy --quiet $(shell find . -type f -name '*.[ch]')

format:
	@clang-format -i $(shell find . -type f -name '*.[ch]')

# dependencies
CFLAGS += -MMD -MP
DEP = $(OBJ:.o=.d) $(BIN:=.d)
-include $(DEP)

# configuration stamp
STAMP := obj/.config-$(CONFIG)
$(STAMP):
	@mkdir -p $(dir $@)
	@rm -f obj/.config-*
	@touch $@

# suffix rules
.SUFFIXES:
obj/%.o: src/%.c $(STAMP) Makefile
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) -c $< -o $@

bin/%: run/%.c $(OBJ)
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) -Wno-unused-parameter $< $(OBJ) $(LDLIBS) -o $@
