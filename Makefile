# compiler and default flags
CC = clang
MPICC = OMPI_MPICC=$(CC) mpicc
CFLAGS = -std=c23 -g -Isrc -Wall -Wextra -Wpedantic -Wshadow

# debug flags
CFLAGS += -Og -fno-omit-frame-pointer

# release flags
#CFLAGS += -DNDEBUG -march=native -Ofast -flto=auto

# profiling flags
#CFLAGS += -pg

# libraries
LDLIBS = -lm -lgmsh -lmetis -lhdf5

# sources, objects, and programs
SRC = $(shell find src -type f -name '*.c')
RUN = $(shell find run -type f -name '*.c')
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))
BIN = $(patsubst run/%.c, bin/%, $(RUN))

# make functions
.PHONY: all clean check format tidy
all: $(OBJ) $(BIN)

clean:
	@rm -rf obj bin

check:
	@cppcheck --quiet --project=compile_commands.json \
		--enable=all --inconclusive --check-level=exhaustive \
		--suppress=checkersReport --suppress=missingIncludeSystem --suppress=unusedFunction

format:
	@clang-format -i $(shell find . -type f -name '*.[ch]')

tidy:
	@clang-tidy --quiet $(shell find . -type f -name '*.[ch]')

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
	-@$(MPICC) $(CFLAGS) $< $(OBJ) $(LDLIBS) -o $@
