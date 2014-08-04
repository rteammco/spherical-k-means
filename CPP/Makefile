# Makefile for Spherical K-Means code


# Galois path (needs to be installed and specified here)
GALOIS_PATH = ./../../../Galois-2.2.1


# specify flags and compiler variables
EXE = spkmeans
FLAGS = -DNDEBUG -g -O3 -ffast-math -fopenmp -std=c++0x
SRC_DIR = src
OBJ_DIR = obj


# specify source files
SRC_FILES = main.cpp reader.cpp vectors.cpp timer.cpp cluster_data.cpp spkmeans.cpp spkmeans_openmp.cpp
GALOIS_SRC_FILES = spkmeans_galois.cpp
OBJ = $(addprefix obj/, $(SRC_FILES:.cpp=.o))
GALOIS_OBJ = $(addprefix obj/, $(GALOIS_SRC_FILES:.cpp=.o))
SRC = $(addprefix $(SRC_DIR)/, $(SRC_FILES))
GALOIS_SRC = $(addprefix $(SRC_DIR)/, $(GALOIS_SRC_FILES))


# specify libraries and include sources
GALOIS_INCLUDES = -I$(GALOIS_PATH)/build/release/include -I$(GALOIS_PATH)/include
GALOIS_LINKS = -L$(GALOIS_PATH)/build/release/src -lgalois
LINKS = -pthread -ldl


# set a flag if the Galois library path exists and if galois was installed
GALOIS_EXISTS  = $(shell if [ -d $(GALOIS_PATH) ]; then echo y; fi)
GALOIS_EXISTS += $(shell if [ -d $(GALOIS_PATH)/build/release/lib ]; then echo y; fi)


################################################################################
################################################################################


all: prep compile

prep:
	@mkdir -p $(OBJ_DIR)

# if Galois is installed, compile everything:
ifeq ($(GALOIS_EXISTS), y y)
$(OBJ): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $@..."
	@g++ $(FLAGS) $(GALOIS_INCLUDES) -c -o $@ $<
$(GALOIS_OBJ): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $@..."
	@g++ $(FLAGS) $(GALOIS_INCLUDES) -c -o $@ $<

compile: $(OBJ) $(GALOIS_OBJ)
	@echo "\nBuilding executable:"
	g++ $(FLAGS) $(GALOIS_INCLUDES) $? -o $(EXE) $(GALOIS_LINKS) $(LINKS)
	@echo "\nCompilation complete: success!\n"

# otherwise, if Galois is not installed, compile without it:
# TODO - else ifeq($(GALIS_EXISTS), y) ... only works in new version of make
else
	@echo "\nGalois library not found... compiling without it.\n"
$(OBJ): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $@..."
	@g++ $(FLAGS) -c -o $@ $< -D NO_GALOIS

compile: $(OBJ)
	@echo "\nBuilding executable:"
	g++ $(FLAGS) $? -o $(EXE) $(LINKS) -D NO_GALOIS
	@echo "\nCompilation complete: success!\n"
endif

################################################################################

# Remove the executable and object files
clean:
	@echo "Deleting object files and executable:"
	rm -f spkmeans
	rm -rf $(OBJ_DIR)


################################################################################
################################################################################


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
