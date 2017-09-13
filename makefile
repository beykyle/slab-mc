OBJS = Transport.o Random.o
CC = g++
DEBUG = -g
CFLAGS = -O1 -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

runTransport : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o runTransport

Transport.o : Transport.cpp Random.h
	$(CC) $(CFLAGS) Transport.cpp 

RANDOM.O :  Random.cpp
	$(CC) $(CFLAGS) Random.cpp
clean:
	rm *.o *~ runTransport

tar:
	tar cfv runTransport.tar Transport.cpp Random.cpp Random.h

