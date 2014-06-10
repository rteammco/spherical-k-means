# Include the Galois library (relative path):
GALOIS_PATH = ./../../../Galois-2.2.1

INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois -pthread -ldl 
SRC = src/main.cpp src/reader.cpp src/vectors.cpp src/spkmeans.cpp src/spkmeans_openmp.cpp src/spkmeans_galois.cpp


# Compile normally (includes debug information):
all:
	g++ -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x $(SRC) -o spkmeans $(INCLUDES) $(LINKS)


# Remove the executable:
clean:
	rm spkmeans


# Run with a small test file:
test:
	./spkmeans -d test.txt


# Run with the larger test data set:
run:
	./spkmeans -d ../TestData/documents -v ../TestData/vocabulary


# Debug with GDB (running on larger test data set):
#debug:
#	gdb ./spkmeans ../TestData/documents
