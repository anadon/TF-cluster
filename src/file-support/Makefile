OBJECTS=file-support.o

all: $(OBJECTS)

$(OBJECTS):
	cd TF-Cluster ; make
	ld -r TF-Cluster/*.o -o $(OBJECTS)

clean:
	rm -f $(FILESUPPORT)
	cd TF-Cluster ; make clean
