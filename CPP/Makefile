# Include the Galois library (relative path):
GALOIS_PATH = ./../../../Galois-2.2.1

# specify libraries and include files
GALOIS_INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
GALOIS_LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois
LINKS = -pthread -ldl
SRC = src/main.cpp src/reader.cpp src/vectors.cpp src/timer.cpp src/cluster_data.cpp src/spkmeans.cpp src/spkmeans_openmp.cpp
GALOIS_SRC = src/spkmeans_galois.cpp


# set a flag if the Galois library path exists and if galois was installed
GALOIS_EXISTS =  $(shell if [ -d $(GALOIS_PATH) ]; then echo y; fi)
GALOIS_EXISTS += $(shell if [ -d $(GALOIS_PATH)/build/release/lib ]; then echo y; fi)


# compile commands: with and without Galois
COMPILE = g++ -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x $(SRC) $(GALOIS_SRC) -o spkmeans $(GALOIS_INCLUDES) $(GALOIS_LINKS) $(LINKS)
COMPILE_NO_GALOIS = g++ -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x $(SRC) -o spkmeans $(LINKS) -D NO_GALOIS


# Compile normally (includes debug information):
all:
ifeq ($(GALOIS_EXISTS), y y)
	$(COMPILE)
# TODO - this only works with new versions of make
#else ifeq($(GALOIS_EXISTS), y)
#	@echo "Galois library found but not installed. Compiling without Galois."
#	$(COMPILE_NO_GALOIS)
else
	@echo "Galois library not found. Compiling without Galois."
	$(COMPILE_NO_GALOIS)
endif


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
COMMAND4 = -d ../TestData/news20 -k 200 --noresults
comp4:
	./spkmeans $(COMMAND4); ./spkmeans2 $(COMMAND4)


# Debug with GDB (running on larger test data set):
#debug:
#	gdb ./spkmeans ../TestData/documents
