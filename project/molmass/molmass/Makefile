MolarMass.o: MolarMass.cpp MolarMass.h
	g++ -c -O2 -I. -o MolarMass.o MolarMass.cpp

ListOfMass.o: ListOfMass.cpp ListOfMass.h
	g++ -c -O2 -I. -o ListOfMass.o ListOfMass.cpp

main.o: main.cpp MolarMass.o ListOfMass.o
	g++ -c -O2 -I. -o main.o main.cpp

all: Makefile MolarMass.o ListOfMass.o main.o
	g++ -o prova_mol main.o MolarMass.o ListOfMass.o

clean:
	rm -f prova_mol main.o MolarMass.o ListOfMass.o
