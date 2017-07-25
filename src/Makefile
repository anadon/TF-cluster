CPP=g++
INCLUDES=../include
FILESUPPORT=file-support/file-support.o
CCOMMON=-pipe -Wall -Wextra -Wconversion -std=c++11 -march=native -I$(INCLUDES) -Wmissing-declarations
CFLAGS=-O3 $(CCOMMON)
LIBS=-pthread -lmadlib

EXEC=tf-cluster


SOURCES=main.cpp auxillaryUtilities.cpp tripleLink.cpp geneData.cpp    \
        diagnostics.cpp
OBJECTS=main.o auxillaryUtilities.o tripleLink.o geneData.o            \
        diagnostics.o

all:$(EXEC)


$(EXEC): $(OBJECTS) $(FILESUPPORT)
	$(CPP) $(CFLAGS) -flto $(OBJECTS) $(FILESUPPORT) $(LIBS)  -o $(EXEC)

%.o:%.cpp $(HEADERS)
	$(CPP) $(CFLAGS) -c $<

$(FILESUPPORT):
	cd file-support; make
	cp -f $(FILESUPPORT) .

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
	rm -f gmon.out
	cd file-support ; make clean
