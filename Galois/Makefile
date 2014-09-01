# Makefile for Galois Spherical K-Means code


# Galois path (needs to be installed and specified here)
GALOIS_PATH = ./../../../Galois-2.2.1


# specify flags and compiler variables
EXE = spkmeans
COMPILER = /usr/bin/g++
FLAGS = -DNDEBUG -g -O3 -ffast-math -std=c++0x
SRC_DIR = src
OBJ_DIR = obj


# specify source files
SRC_FILES = spkmeans_galois.cpp
SRC = $(addprefix $(SRC_DIR)/, $(SRC_FILES))


# specify libraries and include sources
INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois -pthread -ldl


################################################################################
################################################################################


# compile everything
all:
	@echo "\nBuilding executable:"
	$(COMPILER) $(FLAGS) $(INCLUDES) $(SRC) -o $(EXE) $(LINKS)
	@echo "\nCompilation complete: success!\n"


################################################################################


# Remove the executable and object files
clean:
	@echo "Deleting executable: $(EXE)"
	rm -f spkmeans


################################################################################
################################################################################

run:
	./spkmeans ../TestData/mat.gr 2 2
