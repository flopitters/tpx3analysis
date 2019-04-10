# ---------------------------------------------------------------------------
# File: Makefile
# ------------------------
#
# Description:
# Makefile for tpx3Analysis.
#
# Version: 0.7
# Author: Florian Pitters
# ---------------------------------------------------------------------------



# Preparing
# -------------------------------------

# Compiler
CC = clang++

# Directories
BINDIR = bin/
COREDIR = core/
ALGDIR = algorithms/
INCLUDE_PATH = 	-I core \
		-I algorithms

# Compiler flags
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs)
CFLAGS = -Wall ${ROOTCFLAGS} ${INCLUDE_PATH}
LFLAGS = -O3 ${ROOTLIBS}

# Executables
EXE = run

# Automatically decide what to compile
CORE = $(notdir $(wildcard core/*.cc))
OBJS = $(CORE:.cc=.o)
OBJS := $(addprefix ${COREDIR}, ${OBJS})

ALGORITHMS = $(notdir $(wildcard algorithms/*.cc))
ALGOBJS = $(ALGORITHMS:.cc=.o)
ALGOBJS := $(addprefix ${ALGDIR}, ${ALGOBJS})



# Processing
# -------------------------------------

# Compile core, user algorithms and make an executable for each user algorithm
all: ${OBJS} ${ALGOBJS} ${BINDIR}${EXE}
	@echo "Done"

${BINDIR}${EXE}: ${OBJS} ${ALGOBJS}
	@echo "Making executable $(notdir $@)"
	@${CC} ${CFLAGS} ${OBJS} ${ALGOBJS} ${LFLAGS} -o $@

${COREDIR}${EXE}.o: core/${EXE}.cc
	@echo "Compiling $(notdir $<)"
	@${CC} $(CFLAGS) -c $< -o $@

${COREDIR}%.o: core/%.cc core/%.h
	@echo "Compiling $(notdir $<)"
	@${CC} $(CFLAGS) -c $< -o $@

${ALGDIR}%.o: algorithms/%.cc algorithms/%.h
	@echo "Compiling $(notdir $<)"
	@$(CC) $(CFLAGS) -c $< -o $@



# Testing and Cleaning
# -------------------------------------

# Test the executable
test:
	@./bin/run -a hello_world
	@echo "All tests successfull"

# Remove all executables and object files
clean:
	@rm -f ${BINDIR}*
	@rm -f ${COREDIR}*.o
	@rm -f ${ALGDIR}*.o
	@echo "Cleaning"
