# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
DEBUG_FLAGS = -g -DDEBUG -O0 -Wall

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

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
EXECUTABLE = $(BUILD_DIR)/mc_sim

# Default target
all: $(EXECUTABLE)

# Main executable
$(EXECUTABLE): $(ALL_OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Object files compilation
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -I$(SRC_DIR) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(EXECUTABLE)

# Run the program
run: $(EXECUTABLE)
	./$(EXECUTABLE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Phony targets
.PHONY: all debug run clean

# Help target
help:
	@echo "Available targets:"
	@echo "  all    - Build executable (default)"
	@echo "  debug  - Build with debug flags"
	@echo "  run    - Build and run program"
	@echo "  clean  - Remove build directory"
	@echo "  help   - Show this help message"