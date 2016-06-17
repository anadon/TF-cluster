CPP=g++
##CFLAGS=-ggdb -pg -O0 -pipe -Wall -std=c++0x
CFLAGS=-O3 -pipe -flto -Wall -Wextra -std=c++0x -march=native
LIBS=-pthread

EXEC=correlation-matrix-pthread

SOURCES=main.cpp auxillaryUtilities.cpp statistics.cpp
OBJECTS=main.o   auxillaryUtilities.o   statistics.o
HEADERS=         auxillaryUtilities.hpp statistics.hpp

all:$(EXEC)

$(EXEC):$(OBJECTS)
	$(CPP) $(CFLAGS) $(LIBS) $(OBJECTS) -o $(EXEC)

##%.o:%.c $(HEADERS)
##	$(CPP) $(CFLAGS) -c  $< -o $@

main.o:main.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c main.cpp -o main.o

auxillaryUtilities.o:auxillaryUtilities.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c auxillaryUtilities.cpp -o auxillaryUtilities.o

statistics.o:statistics.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c statistics.cpp -o statistics.o

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)