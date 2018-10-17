OBJS 	= k_medoidsFINAL.o Functions.o 
SOURCE	= k_medoidsFINAL.c Functions.c
HEADER  = Functions.h 
OUT  	= medoids
CC	= gcc
CFLAGS   = -g -c 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

k_medoidsFINAL.o: k_medoidsFINAL.c
	$(CC) $(CFLAGS) k_medoidsFINAL.c

Functions.o: Functions.c
	$(CC) $(CFLAGS) Functions.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)