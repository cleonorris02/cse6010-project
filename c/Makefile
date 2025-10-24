CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -pedantic -O2

OBJS=error_detection.o main.o

all: error_detection

error_detection: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

error_detection.o: error_detection.c error_detection.h

main.o: main.c error_detection.h

clean:
	rm -f $(OBJS) error_detection

.PHONY: all clean
