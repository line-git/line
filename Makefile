USE_OPENMP ?= yes

# detect the operating system
ifeq ($(shell uname),Darwin)  # for macOS
	CC = clang
	CXX = c++
	CXXFLAGS = -std=c++17

else ifeq ($(shell uname),Linux)  # for Linux
	CC = gcc
	CXX = g++
	CXXFLAGS = -std=c++17

else  # unsupported operating system
	$(error Unsupported operating system)
endif

ifeq ($(USE_OPENMP), yes)
	CFLAGS += -fopenmp
	CXXFLAGS += -fopenmp
	CXXFLAGS += -DUSE_OPENMP
endif

# CFLAGS += -Wall
# CXXFLAGS += -Wall

# Add debug flags if DEBUG variable is passed via the command line
ifeq ($(DEBUG),-g)
	CFLAGS += -g
	CXXFLAGS += -g
endif

# directories
INCLUDE_DIR = include
SRC_DIR = src
OBJ_DIR = obj

# output file
OUTPUT = line

# source files
C_SOURCES = $(wildcard $(SRC_DIR)/*.c)
CPP_SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

# object files
C_OBJECTS = $(C_SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
CPP_OBJECTS = $(CPP_SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# header files
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)

# libraries
LDFLAGS = -lm -lmpc -lmpfr -lgmp -lcrypto -lssl

# create obj directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

# main rule
$(OUTPUT): $(C_OBJECTS) $(CPP_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# rule for .c files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# rule for .cpp files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# clean generated files
clean:
	rm -rf $(OBJ_DIR) $(OUTPUT)

# rebuild rule
rebuild: clean $(OUTPUT)
