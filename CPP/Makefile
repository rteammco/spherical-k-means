# Include the Galois library (relative path):
GALOIS_PATH = ./../../../Galois-2.2.1

INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois -pthread -ldl 
SRC = src/main.cpp src/reader.cpp src/vectors.cpp src/cluster_data.cpp src/spkmeans.cpp src/spkmeans_openmp.cpp src/spkmeans_galois.cpp


# Compile normally (includes debug information):
all:
	g++ -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x $(SRC) -o spkmeans $(INCLUDES) $(LINKS)


# Remove the executable:
clean:
	rm spkmeans


# Run with a small test file:
test:
	./spkmeans -d ../TestData/documents --autok --noresults


# Run with the larger test data set:
run:
	./spkmeans -d ../TestData/documents --auto --noresults --galois


# Compare two versions
COMMAND = -d ../TestData/documents --auto --noresults
comp:
	./spkmeans $(COMMAND); ./spkmeans2 $(COMMAND)


# Debug with GDB (running on larger test data set):
#debug:
#	gdb ./spkmeans ../TestData/documents
