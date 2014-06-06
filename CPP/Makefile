# Include the Galois library (relative path):
GALOIS_PATH = ./../../../Galois-2.2.1

INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois -pthread -ldl 
SRC = src/spkmeans.cpp


# Compile normally (includes debug information):
all:
	g++ -DNDEBUG -g -O3 -ffast-math -std=c++0x $(SRC) -o spkmeans $(INCLUDES) $(LINKS)


# Remove the executable:
clean:
	rm spkmeans


# Run with a small test file:
test:
	./spkmeans test.txt


# Run with the larger test data set:
run:
	./spkmeans ../TestData/documents


# Debug with GDB (running on larger test data set):
#debug:
#	gdb ./spkmeans ../TestData/documents
