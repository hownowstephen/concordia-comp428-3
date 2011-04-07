CC=mpiCC
LIBS= -lmpi

main: main.cc
	$(CC) $(COPT) $(CFLAGS) main.cc -o parallel $(LIBS)

clean:
	rm parallel
	rm -f core*