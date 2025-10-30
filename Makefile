# compiler, default flags, and libraries
CC = clang
CFLAGS = -std=c99 -g -Wall -Wextra -Wpedantic -Wshadow -Wwrite-strings -Isrc
LDLIBS = -lm -lhdf5 -lmetis -lparmetis -lunwind -lunwind-x86_64
MPICC = OMPI_CC=$(CC) MPICH_CC=$(CC) mpicc

# debug flags
CFLAGS += -O0 -fno-omit-frame-pointer -fsanitize-trap -fsanitize=undefined
CFLAGS += -fsanitize=address
#CFLAGS += -D ENABLE_VALGRIND

# release flags
#CFLAGS += -O3 -march=native -flto=auto -D NDEBUG

# profiling flags
#CFLAGS += -pg -fno-inline-functions -fno-optimize-sibling-calls

# sources, objects, and programs
SRC = $(shell find src -type f -name '*.c')
RUN = $(shell find run -type f -name '*.c')
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))
BIN = $(patsubst run/%.c, bin/%, $(RUN))

# make functions
.PHONY: all clean check tidy format
all: $(BIN)

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

# build rules
.SUFFIXES:
obj/%.o: src/%.c Makefile
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) -c $< -o $@

bin/%: run/%.c $(OBJ)
	@mkdir -p $(@D)
	@$(MPICC) $(CFLAGS) $< $(OBJ) $(LDLIBS) -o $@
