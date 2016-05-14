CPP=g++
CFLAGS=-g -O0 -pipe -Wall -std=c++0x
##LIBS=

EXEC=TF-cluster

SOURCES=main.cpp auxillaryUtilities.cpp nodeNetwork.cpp gEdge.cpp gNode.cpp
OBJECTS=main.o auxillaryUtilities.o nodeNetwork.o gEdge.o gNode.o

all:$(OBJECTS)
	$(CPP) $(CFLAGS) $(LIBS) $(OBJECTS) -o $(EXEC)

main.o:main.cpp pairTable.hpp auxillaryUtilities.hpp
	$(CPP) $(CFLAGS) -c main.cpp

auxillaryUtilities.o:auxillaryUtilities.cpp auxillaryUtilities.hpp

nodeNetwork.o:nodeNetwork.cpp nodeNetwork.hpp auxillaryUtilities.hpp gEdge.hpp gNode.hpp
	$(CPP) $(CFLAGS) -c nodeNetwork.cpp

gEdge.o:gEdge.cpp gEdge.hpp gNode.hpp
	$(CPP) $(CFLAGS) -c gEdge.cpp

gNode.o:gNode.cpp gEdge.hpp gNode.hpp nodeNetwork.hpp
	$(CPP) $(CFLAGS) -c gNode.cpp

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)