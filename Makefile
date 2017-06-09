CPP=g++
#CFLAGS=-ggdb -O0 -pipe -Wall -Wextra -Wconversion -std=c++11 -march=native
CFLAGS=-O3 -pipe -Wall -Wextra -Wconversion -std=c++0x -march=native
LIBS=-pthread
CMTX=correlation-matrix.a

EXEC=tf-cluster

SOURCES=main.cpp auxillaryUtilities.cpp tripleLink.cpp geneData.cpp diagnostics.cpp
OBJECTS=main.o   auxillaryUtilities.o   tripleLink.o   geneData.o diagnostics.o
HEADERS=auxillaryUtilities.hpp edge.hpp geneData.hpp graph.hpp \
        tripleLink.hpp vertex.hpp diagnostics.hpp
CMTX_INCLUDE=correlation-matrix.hpp statistics.h
TEMPLATES=edge.t.hpp graph.t.hpp vertex.t.hpp \
					upper-diagonal-square-matrix.t.hpp

all:$(EXEC)


$(EXEC):$(CMTX) $(OBJECTS)
	$(CPP) $(CFLAGS) -flto $(OBJECTS) $(LIBS) $(CMTX) -o $(EXEC)

%.o:%.cpp $(HEADERS) $(TEMPLATES) $(CMTX_INCLUDE)
	$(CPP) $(CFLAGS) -c $<

$(CMTX_INCLUDE):$(CMTX)
	cp correlation-matrix/correlation-matrix.hpp .
	cp correlation-matrix/statistics.h .

$(CMTX):
	cd correlation-matrix ; make
	cp -f correlation-matrix/correlation-matrix.a .
	cp -f correlation-matrix/correlation-matrix.o .
	cp -f correlation-matrix/statistics.o .

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
	rm -f $(CMTX) $(CMTX_INCLUDE)
	rm -f gmon.out
	cd correlation-matrix/ ; make clean

install:$(EXEC)
	install -C -t /usr/bin/ $(EXEC)
