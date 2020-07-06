CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o individuo.o 
	$(CC) random.o individuo.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
individuo.o : individuo.cpp individuo.h
	$(CC) -c individuo.cpp -o individuo.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
