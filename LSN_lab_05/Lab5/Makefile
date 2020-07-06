CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o position.o
	$(CC) random.o position.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
position.o : position.cpp position.h
	$(CC) -c position.cpp -o position.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
