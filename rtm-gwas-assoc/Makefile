# make
# make win32
# make win64
# make macos

TARGET   := rtm-gwas-assoc

CC       := gcc
CXX      := g++
CFLAGS   := -O2 -std=c11 -Wall -DNDEBUG
CXXFLAGS := -O2 -std=c++11 -fopenmp -Wall -DNDEBUG
LDFLAGS  := -s -static -fopenmp
LDLIBS   := -llapack -lblas -lgfortran -lquadmath

OBJ := $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))
OBJ += $(patsubst %.c,%.o,$(wildcard src/*.c))

ifeq ($(MAKECMDGOALS),win32)
CC  := i686-w64-mingw32-gcc
CXX := i686-w64-mingw32-g++
EXE := $(addsuffix .exe,$(TARGET))
endif

ifeq ($(MAKECMDGOALS),win64)
CC  := x86_64-w64-mingw32-gcc
CXX := x86_64-w64-mingw32-g++
EXE := $(addsuffix .exe,$(TARGET))
endif

ifeq ($(MAKECMDGOALS),macos)
CC  := clang
CXX := clang++
LDFLAGS := $(filter-out -s -static,$(LDFLAGS))
LDLIBS  := -framework Accelerate
LDFLAGS += -L/usr/local/lib
CXXFLAGS += -I/usr/local/include
endif

EXE := $(TARGET)

.PHONY: all win32 win64 macos clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $? $(LDLIBS)

win32: all

win64: all

macos: all

clean:
	$(RM) $(OBJ) $(TARGET) $(TARGET).exe
