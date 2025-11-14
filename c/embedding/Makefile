/*         
 * Outline and logic generated with ChatGPT (OpenAI), Nov 2025.         
 * Reviewed and modified by Cleo Norris.         
*/

CC=clang
CFLAGS ?= -std=c11 -Wall -Wextra -Werror -pedantic
INCLUDES := -Ipython

LIB := python/libembedding.a
DEMO := python/embedding_demo

all: $(DEMO)

$(DEMO): $(LIB) python/demo_embedding.o
	$(CC) $(CFLAGS) python/demo_embedding.o -Lpython -lembedding -o $@

$(LIB): python/embedding.o
	ar rcs $@ $^

python/embedding.o: python/embedding.c python/embedding.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

python/demo_embedding.o: python/demo_embedding.c python/embedding.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f python/*.o $(LIB) $(DEMO)

.PHONY: all clean
