CPP=g++

INCLUDES=../../../include
CPPFLAGS=-ggdb -pg -O0 -pipe -Wall -Wextra -Wconversion -std=c++0x -march=native
#CPPFLAGS=-O3 -pipe -Wall -Wextra -Wconversion -std=c++0x -march=native

SOURCES= tf-cluster-files.cpp
OBJECTS= tf-cluster-files.o

all: $(OBJECTS)

%.o:%.cpp $(HEADERS)
	$(CPP) $(CPPFLAGS) -I$(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS)
	rm -f gmon.out
