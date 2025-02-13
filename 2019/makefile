SRC_DIRECTORY=src
BIN_DIRECTORY=bin
DOC_DIRECTORY=doc

C_DIRECTORY=C

SMALL_PARTIAL=168
SMALL_PARTIAL_HYBRID=336
SMALL_GLOBAL=672
SMALL_DEFINES=-DROWS=$(SMALL_GLOBAL) -DCOLUMNS=$(SMALL_GLOBAL)
SMALL_DEFINES_MPI_C=-DROWS=$(SMALL_PARTIAL) -DROWS_GLOBAL=$(SMALL_GLOBAL) -DCOLUMNS=$(SMALL_GLOBAL)
SMALL_DEFINES_MPI_FORTRAN=-DROWS=$(SMALL_GLOBAL) -DCOLUMNS_GLOBAL=$(SMALL_GLOBAL) -DCOLUMNS=$(SMALL_PARTIAL)
SMALL_DEFINES_HYBRID_C=-DROWS=$(SMALL_PARTIAL_HYBRID) -DROWS_GLOBAL=$(SMALL_GLOBAL) -DCOLUMNS=$(SMALL_GLOBAL)
SMALL_DEFINES_HYBRID_FORTRAN=-DROWS=$(SMALL_GLOBAL) -DCOLUMNS_GLOBAL=$(SMALL_GLOBAL) -DCOLUMNS=$(SMALL_PARTIAL_HYBRID)

BIG_PARTIAL=130
BIG_PARTIAL_HYBRID=1820
BIG_GLOBAL=14560
BIG_DEFINES=-DROWS=$(BIG_GLOBAL) -DCOLUMNS=$(BIG_GLOBAL)
BIG_DEFINES_MPI_C=-DROWS=$(BIG_PARTIAL) -DROWS_GLOBAL=$(BIG_GLOBAL) -DCOLUMNS=$(BIG_GLOBAL)
BIG_DEFINES_MPI_FORTRAN=-DROWS=$(BIG_GLOBAL) -DCOLUMNS_GLOBAL=$(BIG_GLOBAL) -DCOLUMNS=$(BIG_PARTIAL)
BIG_DEFINES_HYBRID_C=-DROWS=$(BIG_PARTIAL_HYBRID) -DROWS_GLOBAL=$(BIG_GLOBAL) -DCOLUMNS=$(BIG_GLOBAL)
BIG_DEFINES_HYBRID_FORTRAN=-DROWS=$(BIG_GLOBAL) -DCOLUMNS_GLOBAL=$(BIG_GLOBAL) -DCOLUMNS=$(BIG_PARTIAL_HYBRID)

CC=gcc
MPICC=mpicc
CFLAGS=-std=c99 -O3 -lm

default: quick_compile

all: help documentation quick_compile 

quick_compile: help create_directories serial_versions openmp_versions mpi_versions hybrid_cpu_versions clean_objects

################
# SERIAL CODES #
################
serial_versions: print_serial_compilation C_serial_small C_serial_big

print_serial_compilation:
	@echo -e "\n/////////////////////////////"; \
	 echo "// COMPILING SERIAL CODES //"; \
	 echo "///////////////////////////";

C_serial_small: $(SRC_DIRECTORY)/$(C_DIRECTORY)/serial.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Small-grid version ($(SMALL_GLOBAL)x$(SMALL_GLOBAL))\n        \c";
	$(CC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/serial_small $(SRC_DIRECTORY)/$(C_DIRECTORY)/serial.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(SMALL_DEFINES) -DVERSION_RUN=\"serial_small\"

C_serial_big: $(SRC_DIRECTORY)/$(C_DIRECTORY)/serial.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Big-grid version ($(BIG_GLOBAL)x$(BIG_GLOBAL))\n        \c";
	$(CC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/serial_big $(SRC_DIRECTORY)/$(C_DIRECTORY)/serial.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(BIG_DEFINES) -DVERSION_RUN=\"serial_big\"

################
# OPENMP CODES #
################
openmp_versions: print_openmp_compilation C_openmp_small C_openmp_big

print_openmp_compilation:
	@echo -e "\n/////////////////////////////"; \
	 echo "// COMPILING OPENMP CODES //"; \
	 echo "///////////////////////////";

C_openmp_small: $(SRC_DIRECTORY)/$(C_DIRECTORY)/openmp.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Small-grid version ($(SMALL_GLOBAL)x$(SMALL_GLOBAL))\n        \c";
	$(CC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/openmp_small $(SRC_DIRECTORY)/$(C_DIRECTORY)/openmp.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(SMALL_DEFINES) -DVERSION_RUN=\"openmp_small\" -fopenmp

C_openmp_big: $(SRC_DIRECTORY)/$(C_DIRECTORY)/openmp.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Big-grid version ($(BIG_GLOBAL)x$(BIG_GLOBAL))\n        \c";
	$(CC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/openmp_big $(SRC_DIRECTORY)/$(C_DIRECTORY)/openmp.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(BIG_DEFINES) -DVERSION_RUN=\"openmp_big\" -fopenmp

#############
# MPI CODES #
#############
mpi_versions: print_mpi_compilation C_mpi_small C_mpi_big
print_mpi_compilation:
	@echo -e "\n//////////////////////////"; \
	 echo "// COMPILING MPI CODES //"; \
	 echo "////////////////////////";

C_mpi_small: $(SRC_DIRECTORY)/$(C_DIRECTORY)/mpi.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Small-grid version ($(SMALL_GLOBAL)x$(SMALL_GLOBAL))\n        \c";
	$(MPICC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/mpi_small $(SRC_DIRECTORY)/$(C_DIRECTORY)/mpi.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(SMALL_DEFINES_MPI_C) -DVERSION_RUN=\"mpi_small\" -DVERSION_RUN_IS_MPI

C_mpi_big: $(SRC_DIRECTORY)/$(C_DIRECTORY)/mpi.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Big-grid version ($(BIG_GLOBAL)x$(BIG_GLOBAL))\n        \c";
	$(MPICC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/mpi_big $(SRC_DIRECTORY)/$(C_DIRECTORY)/mpi.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(BIG_DEFINES_MPI_C) -DVERSION_RUN=\"mpi_big\" -DVERSION_RUN_IS_MPI

####################
# HYBRID CPU CODES #
####################
hybrid_cpu_versions: print_hybrid_cpu_compilation C_hybrid_cpu_small C_hybrid_cpu_big

print_hybrid_cpu_compilation:
	@echo -e "\n/////////////////////////////////"; \
	 echo "// COMPILING HYBRID CPU CODES //"; \
	 echo "///////////////////////////////";

C_hybrid_cpu_small: $(SRC_DIRECTORY)/$(C_DIRECTORY)/hybrid_cpu.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Small-grid version ($(SMALL_GLOBAL)x$(SMALL_GLOBAL))\n        \c";
	$(MPICC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/hybrid_cpu_small $(SRC_DIRECTORY)/$(C_DIRECTORY)/hybrid_cpu.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(SMALL_DEFINES_HYBRID_C) -fopenmp -DVERSION_RUN=\"hybrid_cpu_small\" -DVERSION_RUN_IS_MPI

C_hybrid_cpu_big: $(SRC_DIRECTORY)/$(C_DIRECTORY)/hybrid_cpu.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c
	@echo -e "    - [C] Big-grid version ($(BIG_GLOBAL)x$(BIG_GLOBAL))\n        \c";
	$(MPICC) -o $(BIN_DIRECTORY)/$(C_DIRECTORY)/hybrid_cpu_big $(SRC_DIRECTORY)/$(C_DIRECTORY)/hybrid_cpu.c $(SRC_DIRECTORY)/$(C_DIRECTORY)/util.c $(CFLAGS) $(BIG_DEFINES_HYBRID_C) -fopenmp -DVERSION_RUN=\"hybrid_cpu_big\" -DVERSION_RUN_IS_MPI


clean_objects:
	@rm -f *.o *.mod;

#############
# UTILITIES #
#############
# verify_modules:
# 	@if ! type "pgcc" > /dev/null 2>&1; then \
# 		clear; \
# 		echo -e "\n    . "; \
# 		echo -e "   / \\"; \
# 		echo -e "  / ! \\  It looks like the PGI compiler is not loaded."; \
# 		echo -e " /_____\\ On Bridges please issue 'module load cuda/9.2 mpi/pgi_openmpi/19.4-nongpu'. You can now make again :)\n"; \
# 		exit -1; \
# 	fi

create_directories:
	@if [ ! -d $(BIN_DIRECTORY) ]; then mkdir $(BIN_DIRECTORY); fi; \
	if [ ! -d $(BIN_DIRECTORY)/$(C_DIRECTORY) ]; then mkdir $(BIN_DIRECTORY)/$(C_DIRECTORY); fi; \
	if [ ! -d $(BIN_DIRECTORY)/$(FORTRAN_DIRECTORY) ]; then mkdir $(BIN_DIRECTORY)/$(FORTRAN_DIRECTORY); fi 

help:
	@clear; \
	echo "+-----------+"; \
	echo "| Quick help \\"; \
	echo "+-------------+------------------+--------------------+"; \
	echo "| Generate the documentation     | make documentation |"; \
	echo "| Delete all binaries            | make clean         |"; \
	echo -e "+-----------------------------------------------------+\n";

clean: help
	@echo -e "\n////////////////////////";
	@echo "// CLEANING BINARIES //";
	@echo "//////////////////////";
	rm -rf $(BIN_DIRECTORY);

documentation: help
	@echo -e "\n///////////////////////////////";
	@echo "// GENERATING DOCUMENTATION //";
	@echo "/////////////////////////////";
	@echo -e "    - Generating doxygen... \c"; \
	 doxygen > /dev/null 2>&1; \
	 echo "done"; \
	 echo "    - The HTML documentation is available in 'doc/html/index.xhtml'."; \
