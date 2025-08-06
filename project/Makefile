# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -DNDEBUG -O3 -Wall -fopenmp
DEBUG_FLAGS = $(CXXFLAGS) -g -DDEBUG
PROFILER_FLAGS = $(CXXFLAGS) -pg -g
GPROF_FLAGS = --demangle --display-unused-functions

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
DATA_DIR = data/config
PROFILER_DIR = profiler

# Source files
CORE_SOURCES = $(wildcard $(SRC_DIR)/core/*.cpp)
UTILS_SOURCES = $(wildcard $(SRC_DIR)/utils/*.cpp)
MAIN_SOURCES = $(wildcard $(SRC_DIR)/main/*.cpp)

# Object files
CORE_OBJECTS = $(CORE_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
UTILS_OBJECTS = $(UTILS_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
MAIN_OBJECTS = $(MAIN_SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# All object files
ALL_OBJECTS = $(CORE_OBJECTS) $(UTILS_OBJECTS) $(MAIN_OBJECTS)

# Executable
EXEC = $(BUILD_DIR)/mc_sim

# Json file
DATA_FILE = $(DATA_DIR)/simulation.json

# Main executable
$(EXEC): $(ALL_OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Object files compilation
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -I$(SRC_DIR) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Phony targets
.PHONY: build debug run clean profiler help

# Build target
all: $(EXEC)	

# Run the program
run: 
	./$(EXEC) $(DATA_FILE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Debug build
debug:
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory all CXXFLAGS="$(DEBUG_FLAGS)"
	@$(MAKE) --no-print-directory run
 
# Profiler build and run
profile: 
	rm -rf $(PROFILER_DIR)
	mkdir -p $(PROFILER_DIR)
	@$(MAKE) --no-print-directory clean
	@$(MAKE) --no-print-directory all CXXFLAGS="$(PROFILER_FLAGS)"
	@$(MAKE) --no-print-directory run
	gprof $(GPROF_FLAGS) $(EXEC) > $(PROFILER_DIR)/report.txt
	rm -f gmon.out

# Help target
help:
	@echo "Available targets:"
	@echo "  all  		 	- Build executable (default)"
	@echo "  run     		- Run program using data from 'data/config/simulation.json'"
	@echo "  clean  		- Remove build directory"
	@echo "  debug  		- Build and run with debug flags"
	@echo "  profile 		- Build and run with profiling enabled, gprof report saved in 'profiler/'"
	@echo "  help   		- Show this help message"