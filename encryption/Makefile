CC ?= gcc
CFLAGS ?= -O2 -Wall -Wextra -std=c11
LDFLAGS ?=
LIBSODIUM_CFLAGS ?=
LIBSODIUM_LDFLAGS ?= -lsodium
OPENMP_FLAG ?= -fopenmp

SRC_DIR := src
BIN_DIR := bin
INC_DIR := include

SOURCES := $(wildcard $(SRC_DIR)/*.c)
OBJECTS := $(patsubst $(SRC_DIR)/%.c,$(BIN_DIR)/%.o,$(SOURCES))
TARGET := $(BIN_DIR)/hotspot_encrypt

.PHONY: all clean

all: $(TARGET)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/%.o: $(SRC_DIR)/%.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $(LIBSODIUM_CFLAGS) -I$(INC_DIR) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OPENMP_FLAG) $^ -o $@ $(LIBSODIUM_LDFLAGS) $(LDFLAGS)

clean:
	rm -rf $(BIN_DIR)
