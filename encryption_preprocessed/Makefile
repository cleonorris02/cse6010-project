# Compiler (Apple Clang by default)
CC ?= clang

# Homebrew prefixes
BREW_PREFIX       := $(shell brew --prefix)
LIBOMP_PREFIX     := $(shell brew --prefix libomp)
SODIUM_PREFIX     := $(shell brew --prefix libsodium)

# Compiler and linker flags
CFLAGS  ?= -O2 -Wall -Wextra -std=c11 -Xpreprocessor -fopenmp
INCLUDES = -Iinclude -I$(LIBOMP_PREFIX)/include -I$(SODIUM_PREFIX)/include
LDFLAGS ?= -L$(LIBOMP_PREFIX)/lib -L$(SODIUM_PREFIX)/lib
LIBS    ?= -lsodium -lomp

# Sources and targets
SOURCES = src/main.c src/sequence.c
OBJECTS = $(SOURCES:.c=.o)
TARGET  = dna_hotspot_encryptor

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
