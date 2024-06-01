# compiler and default flags
CC = clang
MPICC = OMPI_MPICC=$(CC) mpicc
CFLAGS = -std=c2x -g -Isrc

# warning flags
CFLAGS += -Wall -Wextra -Wpedantic -Wshadow -Wfloat-equal -Wcast-qual
CFLAGS += -Wno-gnu-zero-variadic-macro-arguments

# debug flags
CFLAGS += -Og -fno-omit-frame-pointer

# release flags
#CFLAGS += -march=native -Ofast -flto=auto

# profiling flags
#CFLAGS += -pg -fno-lto -fno-inline

# libraries
LDLIBS = -lm -lmetis -lhdf5 -lgmsh

# sources, objects, and programs
SRC = $(shell find src -type f -name "*.c")
RUN = $(shell find run -type f -name "*.c")
OBJ = $(SRC:src/%.c=obj/%.o)
BIN = $(RUN:run/%.c=bin/%)

# make functions
.PHONY: all clean
all: $(OBJ) $(BIN)

clean:
	rm -rf obj bin

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
