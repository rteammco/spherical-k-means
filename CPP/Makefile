# Include the Galois library (relative path):
GALOIS_PATH = ./../../../Galois-2.2.1

INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois -pthread -ldl 
SRC = src/main.cpp src/reader.cpp src/vectors.cpp src/timer.cpp src/cluster_data.cpp src/spkmeans.cpp src/spkmeans_openmp.cpp src/spkmeans_galois.cpp


# Compile normally (includes debug information):
all:
	g++ -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x $(SRC) -o spkmeans $(INCLUDES) $(LINKS)


# Remove the executable:
clean:
	rm spkmeans


# Run test with a small test file, with different versions:
TESTCMD = ./spkmeans -d ../TestData/documents --autok --noresults
test:
	$(TESTCMD)
# OpenMP version:
testo:
	$(TESTCMD) --openmp
# Galois version:
testg:
	$(TESTCMD) --galois


# Run with the larger test data set:
run:
	./spkmeans -d ../TestData/documents --auto --noresults --galois


# Compare two versions with small dataset
COMMAND = -d ../TestData/documents --auto --noresults
comp:
	./spkmeans $(COMMAND); ./spkmeans2 $(COMMAND)


# Compare two versions w/ big dataset
COMMAND2 = -d ../TestData/classic3 --auto --noresults --noscheme
comp2:
	./spkmeans $(COMMAND2); ./spkmeans2 $(COMMAND2)


# Compare two versions w/ bigger dataset
COMMAND3 = -d ../TestData/mat -k 50 --noresults --noscheme
comp3:
	./spkmeans $(COMMAND3); ./spkmeans2 $(COMMAND3)


# Compare two versions w/ biggest dataset
COMMAND4 = -d ../TestData/news20 -k 10 --noresults
comp4:
	./spkmeans $(COMMAND4); ./spkmeans2 $(COMMAND4)


# Debug with GDB (running on larger test data set):
#debug:
#	gdb ./spkmeans ../TestData/documents
