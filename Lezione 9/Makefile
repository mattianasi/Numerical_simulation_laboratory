CC = g++
CFLAGS = -Wall -O3 --std=c++11
LDFLAGS = -larmadillo

EXEC = main.exe
OBJ = main.o random.o popolazione.o percorso.o

$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(LDFLAGS)

main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

popolazione.o: popolazione.cpp popolazione.h percorso.h random.h
	$(CC) -c popolazione.cpp -o popolazione.o $(CFLAGS)

percorso.o: percorso.cpp percorso.h random.h
	$(CC) -c percorso.cpp -o percorso.o $(CFLAGS)

clean:
	rm -f *.o $(EXEC) seed.out