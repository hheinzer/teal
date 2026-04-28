# compiler
CC = gcc
MPICC = OMPI_CC=$(CC) mpicc

# libraries
LDLIBS = -lm -lmetis -lparmetis -lhdf5

# default flags
CFLAGS = -Isrc -std=c99 -g3 -Wall -Wextra -Wpedantic -Wshadow -Wwrite-strings \
		 -Wno-unused-parameter -Wno-unused-function

# exceptions
ifeq ($(CC), gcc)
	CFLAGS += -Wno-discarded-qualifiers
endif

# debug flags
CFLAGS += -O0 -fno-omit-frame-pointer -fsanitize=address,undefined

# release flags
#CFLAGS += -O3 -march=native -flto=auto -DNDEBUG -Wno-unused

# perf flags
#CFLAGS += -fno-omit-frame-pointer -fno-inline-functions

# sources, objects, and programs
SRC := $(shell find src -type f -name '*.c')
RUN := $(shell find run -type f -name '*.c')
TST := $(shell find test -type f -name '*.c')
OBJ := $(patsubst src/%.c, obj/%.o, $(SRC))
BIN := $(patsubst run/%.c, bin/%, $(RUN)) \
	   $(patsubst test/%.c, bin/test/%, $(TST))

# make functions
.PHONY: all clean check tidy format test

all: $(BIN)

clean:
	@rm -rf obj bin

check:
	@cppcheck -q --project=compile_commands.json --check-level=exhaustive --enable=all \
		--suppress=checkersReport --suppress=missingIncludeSystem \
		--suppress=constVariablePointer --suppress=constVariable --suppress=unusedFunction

tidy: $(OBJ)
	@clang-tidy --quiet $(shell find . -type f -name '*.[ch]')

format:
	@clang-format -i $(shell find . -type f -name '*.[ch]')

test: $(filter bin/test/%, $(BIN))
	@for bin in $(sort $^); do \
		printf "%-30s " "$$bin"; \
		mpirun --oversubscribe -n 8 $$bin -v > "$$bin.out" 2>&1 \
			&& echo "passed" || echo "FAILED"; \
	done

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
