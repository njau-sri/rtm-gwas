TARGET = rtm-gwas-snpldb

RTM_GWAS_VERSION ?= $(shell cat ../VERSION)
DEFINES := -DRTM_GWAS_VERSION=\"$(RTM_GWAS_VERSION)\"

CC  := gcc
CXX := g++

CFLAGS   := -O2 -std=c11 -fopenmp
CXXFLAGS := -O2 -std=c++11 -fopenmp $(DEFINES)
LDFLAGS  := -s -fopenmp

OBJ := $(patsubst %.c,%.o,$(wildcard src/*.c))
OBJ += $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	$(RM) $(OBJ) $(TARGET)
