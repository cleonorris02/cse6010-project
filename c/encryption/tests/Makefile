# Makefile for XChaCha20 verification tests

CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11
LDFLAGS ?=
LIBS ?= -lsodium

# Sources and targets
TEST_SOURCES = test_xchacha20.c
TEST_TARGET = test_xchacha20

.PHONY: all clean test

all: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_SOURCES)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(TEST_TARGET)
