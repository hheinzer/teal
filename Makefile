# compiler and default flags
CC = gcc
MPICC = OMPI_MPICC=$(CC) mpicc
CFLAGS = -std=c2x -g -Isrc

# warning flags
CFLAGS += -Wall -Wextra -Wpedantic -Wshadow -Wfloat-equal -Wcast-qual

# debug flags
CFLAGS += -Og -fno-omit-frame-pointer -fanalyzer

# release flags
#CFLAGS += -march=native -Ofast -flto=auto -DNDEBUG

# profiling flags
#CFLAGS += -pg -fno-lto -fno-inline

# libraries
LDLIBS = -lm -lgmsh -lmetis -lhdf5

# sources, objects, and programs
SRC = $(shell find src -type f -name "*.c")
RUN = $(shell find run -type f -name "*.c")
OBJ = $(SRC:src/%.c=obj/%.o)
BIN = $(RUN:run/%.c=bin/%)

# make functions
.PHONY: all clean check
all: $(OBJ) $(BIN)

clean:
	rm -rf obj bin

check:
	-cppcheck --project=compile_commands.json --enable=all --inconclusive --check-level=exhaustive \
		--suppress=missingIncludeSystem --suppress=unusedFunction

# dependencies
CFLAGS += -MMD -MP
DEP = $(OBJ:.o=.d) $(BIN:=.d)
-include $(DEP)

# build rules
.SUFFIXES:
obj/%.o: src/%.c Makefile
	@mkdir -p $(@D)
	$(MPICC) $(CFLAGS) -c $< -o $@

bin/%: run/%.c $(OBJ) Makefile
	@mkdir -p $(@D)
	-$(MPICC) $(CFLAGS) $< $(OBJ) $(LDLIBS) -o $@
