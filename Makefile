CPP=g++
CFLAGS=-g -O0 -pipe -Wall -std=c++0x
##LIBS=

EXEC=TF-cluster

SOURCES=main.cpp auxillaryUtilities.cpp
OBJECTS=main.o   auxillaryUtilities.o

all:$(OBJECTS)
	$(CPP) $(CFLAGS) $(LIBS) $(OBJECTS) -o $(EXEC)

main.o:main.cpp graph.hpp auxillaryUtilities.hpp
	$(CPP) $(CFLAGS) -c main.cpp

auxillaryUtilities.o:auxillaryUtilities.cpp auxillaryUtilities.hpp
	$(CPP) $(CFLAGS) -c auxillaryUtilities.cpp

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)