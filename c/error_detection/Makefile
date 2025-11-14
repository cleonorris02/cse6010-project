# Outline generated with ChatGPT (OpenAI), Oct 2025.       
# Reviewed and modified by Cleo Norris.

# use Clang compiler
CC=clang

# standard C11, enable all warnings, optimize for speed
CFLAGS=-std=c11 -Wall -Wextra -pedantic -O2

# object files
OBJS=error_detection.o main.o

# default build rule
all: error_detection

# creates an executable file named error_detection
error_detection: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

# tells make how to build error_detection.o
error_detection.o: error_detection.c error_detection.h

# tells make how to build main.o
main.o: main.c error_detection.h

# resets the folder back to source-only state
clean:
	rm -f $(OBJS) error_detection

# mark all and clean as phony targets
.PHONY: all clean