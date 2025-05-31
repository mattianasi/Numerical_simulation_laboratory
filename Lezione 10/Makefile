CXX = mpicxx
CXXFLAGS = -Wall -O3 -std=c++11
LDFLAGS = -larmadillo

EXEC = main.exe
OBJ = main.o random.o popolazione.o percorso.o

$(EXEC): $(OBJ)
	$(CXX) $(OBJ) $(LDFLAGS) -o $(EXEC)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

random.o: random.cpp random.h
	$(CXX) $(CXXFLAGS) -c random.cpp -o random.o

popolazione.o: popolazione.cpp popolazione.h percorso.h random.h
	$(CXX) $(CXXFLAGS) -c popolazione.cpp -o popolazione.o

percorso.o: percorso.cpp percorso.h random.h
	$(CXX) $(CXXFLAGS) -c percorso.cpp -o percorso.o

clean:
	rm -f *.o $(EXEC) seed.out
