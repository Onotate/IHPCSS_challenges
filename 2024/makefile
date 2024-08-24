CC=mpicc
# Set default CFLAGS, allow them to be overridden or extended by the command line
CFLAGS ?= -O2 -Wall -fopenmp -lm
BIN_DIR=bin

# Define source directories
SRC_DIR=src

# Define the directories for binaries
SERIAL_BIN=$(BIN_DIR)/serial
OPENMP_BIN=$(BIN_DIR)/openmp
MPI_BIN=$(BIN_DIR)/mpi
HYBRID_BIN=$(BIN_DIR)/hybrid

default: all

all: serial openmp mpi hybrid

serial: $(SERIAL_BIN)
	$(CC) $(CFLAGS) -o $(SERIAL_BIN)/main_serial $(SRC_DIR)/serial/main.c

openmp: $(OPENMP_BIN)
	$(CC) $(CFLAGS) -fopenmp -o $(OPENMP_BIN)/main_openmp $(SRC_DIR)/openmp/main.c

mpi: $(MPI_BIN)
	$(CC) $(CFLAGS) -o $(MPI_BIN)/main_mpi $(SRC_DIR)/mpi/main.c

hybrid: $(HYBRID_BIN)
	$(CC) $(CFLAGS) -fopenmp -o $(HYBRID_BIN)/main_hybrid $(SRC_DIR)/hybrid/main.c

# Ensure bin directories exist
$(SERIAL_BIN) $(OPENMP_BIN) $(MPI_BIN) $(HYBRID_BIN):
	mkdir -p $@

clean:
	rm -rf $(BIN_DIR)
