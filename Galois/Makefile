# Include the Galois library (relative path):
GALOIS_PATH = ./../../../Galois-2.2.1

INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois -pthread -ldl 
SRC = src/spkmeans.cpp

all:
	g++ -DNDEBUG -g -O3 -ffast-math -std=c++0x $(SRC) -o spkmeans $(INCLUDES) $(LINKS)

clean:
	rm spkmeans

run:
	./spkmeans test.txt
