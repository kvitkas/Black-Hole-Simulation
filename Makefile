CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra

SRC_DIR = src
BUILD_DIR = build
TARGET = black_hole_renderer

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET) realtime_black_hole

realtime: src/main_realtime.cpp
	$(CXX) $(CXXFLAGS) -o realtime_black_hole src/main_realtime.cpp -framework GLUT -framework OpenGL -Wno-deprecated

.PHONY: all clean realtime
