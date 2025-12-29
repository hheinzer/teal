# compiler
CC = gcc
MPICC = OMPI_CC=$(CC) mpicc

# libraries
LDLIBS = -lm -lhdf5 -lmetis -lparmetis

# default flags
CFLAGS = -Isrc -std=c99 -g3 -Wall -Wextra -Wpedantic -Wshadow -Wconversion \
		 -Wno-unused-parameter -Wno-unused-function -Wno-sign-conversion

# debug flags
CFLAGS += -O0 -fno-omit-frame-pointer -fsanitize=address,undefined -fanalyzer

# release flags
#CFLAGS += -O3 -march=native -flto=auto -DNDEBUG

# perf flags
#CFLAGS += -fno-omit-frame-pointer -fno-inline-functions

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
		--suppress=checkersReport --suppress=missingIncludeSystem --suppress=unusedFunction

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
	@$(MPICC) $(CFLAGS) $< $(OBJ) $(LDLIBS) -o $@
