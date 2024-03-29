CC=gcc
CFLAGS=-Wall -std=c99

.c:
	$(CC) $(CFLAGS) -o $@  $<
.o:
	$(CC) $(CFLAGS) -o $@  $<

.c.o:
	$(CC) $(CFLAGS) -c $<

main: arit.o poly.o karatsuba.o toom.o expe.o main.o
	$(CC) $(CFLAGS) main.o -o main

clean:
	rm *.o main