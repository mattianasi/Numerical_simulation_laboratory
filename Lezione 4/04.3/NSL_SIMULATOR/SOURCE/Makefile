CC = g++
CXXFLAGS = -O3 --std=c++14 -I/usr/local/include
LDFLAGS = -L/usr/local/lib -larmadillo

all: simulator.exe

simulator.exe: system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) $(CXXFLAGS) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe $(LDFLAGS)

NSL_SIMULATOR.o: NSL_SIMULATOR.cpp system.h particle.h
	$(CC) $(CXXFLAGS) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o

system.o: system.cpp system.h random.h particle.h
	$(CC) $(CXXFLAGS) -c system.cpp -o system.o

particle.o: particle.cpp particle.h
	$(CC) $(CXXFLAGS) -c particle.cpp -o particle.o

random.o: random.cpp random.h
	$(CC) $(CXXFLAGS) -c random.cpp -o random.o

clean:
	rm -f *.o simulator.exe

remove:
	rm -f ../OUTPUT/*.*
	rm -f ../OUTPUT/CONFIG/*.*
