CC=g++ -march=native -std=c++11 -o3 -fopenmp
CFLAGS=-c -I. -std=c++0x -Wfatal-errors

u_query: u_query.o u_spc.o u_io.o
	$(CC) u_query.o u_spc.o u_io.o -o u_query
	rm *.o

u_index: u_index.o u_spc.o u_io.o
	$(CC) u_index.o u_spc.o u_io.o -o u_index
	rm *.o

u_index.o: u_index.cc
	$(CC) $(CFLAGS) u_index.cc -o u_index.o

u_query.o: u_query.cc
	$(CC) $(CFLAGS) u_query.cc -o u_query.o

u_io.o: u_io.cc
	$(CC) $(CFLAGS) u_io.cc -o u_io.o

u_spc.o: u_spc.cc
	$(CC) $(CFLAGS) u_spc.cc -o u_spc.o
clean:
	rm -rf *.o
