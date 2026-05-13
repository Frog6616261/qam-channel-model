CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -Iinclude -MMD -MP

BUILD_DIR := build
SRC_DIR := src
TARGET := $(BUILD_DIR)/app

SOURCES := $(SRC_DIR)/main.cpp \
           $(SRC_DIR)/awgn.cpp \
           $(SRC_DIR)/qamdemod.cpp \
           $(SRC_DIR)/qammod.cpp

OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DEPS := $(OBJECTS:.o=.d)

.PHONY: all run clean rebuild

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(OBJECTS) -o $(TARGET)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	./$(TARGET)

clean:
	@rm -rf $(BUILD_DIR)

rebuild: clean all

-include $(DEPS)