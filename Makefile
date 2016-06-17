EXEC=correlation-matrix-pthread triple-link-pthread

all:$(EXEC)

correlation-matrix-pthread:
	cd correlation-matrix ; make
	cp correlation-matrix/correlation-matrix-pthread .

triple-link-pthread:
	cd triple-link ; make
	cp triple-link/triple-link-pthread .


clean:
	cd correlation-matrix ; make clean
	cd triple-link ; make clean
	rm -f $(EXEC)