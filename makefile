all: Database.o Sequence.o Algorithm.o
	g++ Database.o Sequence.o Algorithm.o projet.cpp -o projet
Database.o: Database.cpp
	g++ -c Database.cpp
Sequence.o: Database.cpp
	g++ -c Sequence.cpp
Algorithm.o: Sequence.cpp
	g++ -c Algorithm.cpp
clean:
	rm *.o
