CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o funzioni.o
	$(CC) random.o main.o funzioni.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
