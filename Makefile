CPP=g++
##CFLAGS=-g -O0 -pthread -pipe -flto -Wall -Wextra -std=c++0x
CFLAGS=-O3 -pthread -pipe -flto -Wall -Wextra -std=c++0x
##LIBS=

EXEC=TF-cluster

SOURCES=main.cpp auxillaryUtilities.cpp tripleLink.cpp geneData.cpp
OBJECTS=main.o   auxillaryUtilities.o   tripleLink.o   geneData.o

all:$(OBJECTS)
	$(CPP) $(CFLAGS) $(LIBS) $(OBJECTS) -o $(EXEC)

main.o:main.cpp graph.hpp auxillaryUtilities.hpp
	$(CPP) $(CFLAGS) -c main.cpp

auxillaryUtilities.o:auxillaryUtilities.cpp auxillaryUtilities.hpp graph.t.hpp
	$(CPP) $(CFLAGS) -c auxillaryUtilities.cpp

tripleLink.o:tripleLink.cpp tripleLink.hpp graph.hpp
	$(CPP) $(CFLAGS) -c tripleLink.cpp

geneData.o:geneData.cpp geneData.hpp
	$(CPP) $(CFLAGS) -c geneData.cpp

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)