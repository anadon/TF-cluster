CPP=g++
CFLAGS=-ggdb -pg -O0 -pipe -Wall -Wextra -std=c++11 -march=native
##CFLAGS=-O3 -pipe -Wall -Wextra -std=c++0x -march=native
LIBS=-pthread correlation-matrix.a

EXEC=triple-link-pthread

SOURCES=main.cpp auxillaryUtilities.cpp tripleLink.cpp geneData.cpp
OBJECTS=main.o   auxillaryUtilities.o   tripleLink.o   geneData.o
HEADERS=auxillaryUtilities.hpp edge.hpp geneData.hpp graph.hpp \
        tripleLink.hpp vertex.hpp correlation-matrix.hpp
TEMPLATES=edge.t.hpp graph.t.hpp vertex.t.hpp

all:$(EXEC)

$(EXEC):$(OBJECTS)
	$(CPP) $(CFLAGS) -flto $(OBJECTS) $(LIBS) -o $(EXEC)

%.o:%.cpp $(HEADERS) $(TEMPLATES)
	$(CPP) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)