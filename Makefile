#!/bin/make

EXEC=tf-cluster

all:$(EXEC)
	cd src ; make ;  cp $(EXEC) .

clean:
	rm -f $(EXEC)
	rm -f gmon.out
	cd src ; make clean

install:$(EXEC)
	install -C -t /usr/bin/ $(EXEC)
