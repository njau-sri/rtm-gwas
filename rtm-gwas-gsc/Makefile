TARGET = rtm-gwas-gsc

CC  := gcc
CXX := g++

CFLAGS   := -O2 -std=c11
CXXFLAGS := -O2 -std=c++11 -fopenmp
LDFLAGS  := -s -fopenmp
LDLIBS   := -llapack -lblas

OBJ := $(patsubst %.c,%.o,$(wildcard src/*.c))
OBJ += $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $? $(LDLIBS)

.PHONY: clean
clean:
	$(RM) $(OBJ) $(TARGET)
