CPP=g++
#CFLAGS=-ggdb -pg -O0 -pipe -Wall -Wextra -Wconversion -std=c++11 -march=native
CFLAGS=-O3 -pipe -Wall -Wextra -Wconversion -std=c++0x -march=native
LIBS=-pthread
CMTX=correlation-matrix.a

EXEC=triple-link-pthread

SOURCES=main.cpp auxillaryUtilities.cpp tripleLink.cpp geneData.cpp
OBJECTS=main.o   auxillaryUtilities.o   tripleLink.o   geneData.o
HEADERS=auxillaryUtilities.hpp edge.hpp geneData.hpp graph.hpp \
        tripleLink.hpp vertex.hpp
CMTX_INCLUDE=correlation-matrix.hpp statistics.h
TEMPLATES=edge.t.hpp graph.t.hpp vertex.t.hpp

all:$(EXEC)


$(EXEC):$(CMTX) $(OBJECTS)
	$(CPP) $(CFLAGS) -flto $(OBJECTS) $(LIBS) $(CMTX) -o $(EXEC)

%.o:%.cpp $(HEADERS) $(TEMPLATES) $(CMTX_INCLUDE)
	$(CPP) $(CFLAGS) -c $<

$(CMTX_INCLUDE):$(CMTX)
	cp correlation-matrix/correlation-matrix.hpp .
	cp correlation-matrix/statistics.h .

$(CMTX):
	git submodule init
	git submodule update --remote
	cd correlation-matrix ; make
	mv correlation-matrix/correlation-matrix.a .

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
	rm -f $(CMTX) $(CMTX_INCLUDE)
	rm -f gmon.out
	cd correlation-matrix/ ; make clean

install:$(EXEC)
	install -C -t /usr/bin/ $(EXEC)