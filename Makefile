#
# Makefile
#
SOURCES.cpp= Ace.cpp BayesianAlgorithm.cpp BayesianUtils.cpp Contig.cpp DualCodebook.cpp Forage.cpp Key.cpp KeyTable.cpp MeanStdDev.cpp ProbHash.cpp RuntimeParameters.cpp Sequence.cpp SynchronizedArrays.cpp
CC=g++
CPPFLAGS=-g
PROGRAM=forage

OBJECTS= $(SOURCES.cpp:.cpp=.o)

all: $(PROGRAM)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
	$(LINK.cpp) -o $@ $(OBJECTS) $(SLIBS)

clean:
	rm -f $(PROGRAM) $(OBJECTS)
