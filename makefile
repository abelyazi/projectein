all: Database.o
	g++ Database.o projet6.cpp -o projet6
Database.o: Database.cpp
	g++ -c Database.cpp
clean:
	rm *.o
