# Outline and logic generated with ChatGPT (OpenAI), Nov 2025.         
# Reviewed and modified by Cleo Norris.         

CC=clang
CFLAGS ?= -std=c11 -Wall -Wextra -Werror -pedantic
INCLUDES := -Ic

LIB := libembedding.a
DEMO := embedding_demo

all: $(DEMO)

$(DEMO): $(LIB) demo_embedding.o
	$(CC) $(CFLAGS) demo_embedding.o libembedding.a -o $@

$(LIB): embedding.o
	ar rcs $@ $^

c/embedding.o: embedding.c embedding.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

c/demo_embedding.o: demo_embedding.c embedding.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f c/*.o $(LIB) $(DEMO)

.PHONY: all clean
