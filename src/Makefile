# -------- Makefile --------

# ROOT config
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# Compiler and flags
CXX          := g++
CXXFLAGS     := -O2 -Wall -std=c++17

# Sources and outputs
SOURCES      := StParameters.cxx Initiate_input.cpp
OBJECTS      := $(SOURCES:.cxx=.o)
DICT         := StParametersDict.cxx
HEADERS      := StParameters.h
LINKDEF      := StParameters_LinkDef.h
TARGET       := initiate_input_exec

# Default target
all: $(TARGET)

$(DICT): $(HEADERS) $(LINKDEF)
	@echo "Generating dictionary..."
	rootcling -f $@ -c $(HEADERS) $(LINKDEF)

$(TARGET): $(DICT) $(SOURCES)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $(SOURCES) $(DICT) $(ROOTLIBS)

clean:
	rm -f *.o *.pcm *.so *.dylib *.dll $(DICT) $(TARGET)

.PHONY: all clean

