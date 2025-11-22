# compiler
CC = clang
MPICC = OMPI_CC=$(CC) MPICH_CC=$(CC) mpicc

# libraries
LDLIBS = -lm -lhdf5 -lmetis -lparmetis

# default flags
CFLAGS = -Isrc -std=c99 -g -Wall -Wextra -Wpedantic -Wshadow -Wwrite-strings -Wcast-qual

# debug flags
CFLAGS += -O0 -fno-omit-frame-pointer -fsanitize=address,undefined

# valgrind flags
#CFLAGS += -Og -fno-omit-frame-pointer -DVALGRIND

# release flags
#CFLAGS += -O3 -march=native -flto=auto
#CFLAGS += -DNDEBUG -Wno-unused -Wno-unused-parameter

# gprof flags
#CFLAGS += -pg -fno-inline-functions

# sources, objects, and programs
SRC := $(shell find src -type f -name '*.c')
RUN := $(shell find run -type f -name '*.c')
OBJ := $(patsubst src/%.c, obj/%.o, $(SRC))
BIN := $(patsubst run/%.c, bin/%, $(RUN))

# make functions
.PHONY: all clean check tidy format

all: $(BIN)

clean:
	@rm -rf obj bin

check:
	@cppcheck --project=compile_commands.json --check-level=exhaustive --enable=all \
		--suppress=checkersReport --suppress=missingIncludeSystem \
		--suppress=unusedFunction --suppress=constVariable --suppress=constVariablePointer

tidy: $(OBJ)
	@clang-tidy $(shell find . -type f -name '*.[ch]')

format:
	@clang-format -i $(shell find . -type f -name '*.[ch]')

# dependencies
CFLAGS += -MMD -MP
DEP = $(OBJ:.o=.d) $(BIN:=.d)
-include $(DEP)

# suffix rules
.SUFFIXES:
obj/%.o: src/%.c Makefile
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) -c $< -o $@

bin/%: run/%.c $(OBJ)
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) -Wno-unused-parameter $< $(OBJ) $(LDLIBS) -o $@
