# compiler and default flags
CC = clang
MPICC = OMPI_MPICC=$(CC) mpicc
CFLAGS = -std=c23 -g -Isrc

# warning flags
CFLAGS += -Wall -Wextra -Wpedantic -Wshadow -Wfloat-equal -Wcast-qual

# debug flags
CFLAGS += -Og -fno-omit-frame-pointer

# release flags
#CFLAGS += -march=native -Ofast -flto=auto -DNDEBUG

# profiling flags
#CFLAGS += -pg -fno-lto -fno-inline

# libraries
LDLIBS = -lm -lgmsh -lmetis -lhdf5

# sources, objects, and programs
SRC = $(wildcard src/**/*.c)
RUN = $(wildcard run/**/*.c)
OBJ = $(patsubst src/%.c, obj/%.o, $(SRC))
BIN = $(patsubst run/%.c, bin/%, $(RUN))

# make functions
.PHONY: all clean check format tidy
all: $(OBJ) $(BIN)

clean:
	@rm -rf obj bin

check:
	@cppcheck -q --project=compile_commands.json --enable=all --inconclusive --check-level=exhaustive \
		--suppress=checkersReport --suppress=missingIncludeSystem --suppress=unusedFunction

format:
	@clang-format -i $(wildcard src/**/*.[ch]) $(wildcard run/**/*.[ch])

tidy:
	@clang-tidy --quiet $(wildcard src/**/*.[ch]) $(wildcard run/**/*.[ch])

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
