CPP=g++
CFLAGS=-g -O0
LIBS=

EXEC=TF-cluster

SOURCES=main.cpp
OBJECTS=main.o

all:$(OBJECTS)
	$(CPP) $(CFLAGS) $(LIBS) $(OBJECTS) -o $(EXEC)

$(OBJECTS):$(SOURCES)
	$(CPP) $(CFLAGS) -c $<