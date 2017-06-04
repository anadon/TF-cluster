CPP=g++
CFLAGS=-ggdb -O0 -pipe -Wall -Wextra -Wconversion -std=c++11 -march=native
#CFLAGS=-O3 -pipe -Wall -Wextra -Wconversion -std=c++0x -march=native
LIBS=-pthread

EXEC=tf-cluster

SOURCES=main.cpp auxillaryUtilities.cpp tripleLink.cpp geneData.cpp    \
        diagnostics.cpp
OBJECTS=main.o auxillaryUtilities.o tripleLink.o geneData.o            \
        diagnostics.o
HEADERS=auxillaryUtilities.hpp geneData.hpp tripleLink.hpp             \
        diagnostics.hpp

all:$(EXEC)


$(EXEC): $(OBJECTS)
	$(CPP) $(CFLAGS) -flto $(OBJECTS) $(LIBS) -o $(EXEC)

%.o:%.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c $<


clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
	rm -f gmon.out

install:$(EXEC)
	install -C -t /usr/bin/ $(EXEC)
