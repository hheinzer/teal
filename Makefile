# compiler
CC = gcc
MPICC = OMPI_CC=$(CC) mpicc

# libraries
LDLIBS = -lm -lhdf5 -lparmetis

# default flags
CFLAGS = -Isrc -std=c99 -g3 -Wall -Wextra -Wpedantic -Wshadow -Wconversion \
		 -Wno-unused-parameter -Wno-unused-function

# debug flags
CFLAGS += -O0 -fno-omit-frame-pointer -fsanitize=address,undefined -fanalyzer

# release flags
#CFLAGS += -O3 -march=native -flto=auto -DNDEBUG

# perf flags
#CFLAGS += -fno-omit-frame-pointer -fno-inline-functions

# sources, objects, and programs
SRC := $(shell find src -type f -name '*.c')
RUN := $(shell find run -type f -name '*.c')
TEST := $(shell find test -type f -name '*.c')
OBJ := $(patsubst src/%.c, obj/%.o, $(SRC))
BIN := $(patsubst run/%.c, bin/%, $(RUN)) \
       $(patsubst test/%.c, bin/test/%, $(TEST))

# make functions
.PHONY: all check tidy format test clean

all: $(BIN)

check:
	@cppcheck --project=compile_commands.json --check-level=exhaustive --enable=all \
		--suppress=checkersReport --suppress=missingIncludeSystem --suppress=unusedFunction

tidy: $(OBJ)
	@clang-tidy $(shell find . -type f -name '*.[ch]')

format:
	@clang-format -i $(shell find . -type f -name '*.[ch]')

test: $(filter bin/test/%, $(BIN))
	@for exe in $^; do echo $$exe; mpirun -n 2 $$exe -q; done

clean:
	@rm -rf obj bin

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

bin/test/%: test/%.c $(OBJ)
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) $< $(OBJ) $(LDLIBS) -o $@
