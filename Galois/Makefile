# Makefile for Spherical K-Means code


# Galois path (needs to be installed and specified here)
GALOIS_PATH = ./../../../Galois-2.2.1


# specify flags and compiler variables
EXE = spkmeans
COMPILER = /usr/bin/g++
FLAGS = -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x
SRC_DIR = src
OBJ_DIR = obj


# specify source files
GALOIS_SRC_FILES = spkmeans_galois.cpp
GALOIS_SRC = $(addprefix $(SRC_DIR)/, $(GALOIS_SRC_FILES))


# specify libraries and include sources
GALOIS_INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
GALOIS_LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois
LINKS = -pthread -ldl


################################################################################
################################################################################


# compile everything
all:
	@echo "\nBuilding executable:"
	$(COMPILER) $(FLAGS) $(GALOIS_INCLUDES) $(GALOIS_SRC) -o $(EXE) $(GALOIS_LINKS) $(LINKS)
	@echo "\nCompilation complete: success!\n"


################################################################################


# Remove the executable and object files
clean:
	@echo "Deleting executable: $(EXE)"
	rm -f spkmeans


################################################################################
################################################################################

run:
	./spkmeans
