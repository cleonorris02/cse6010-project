## Building the C embedding module

Run `make` from the repository root to compile the reusable library (`c/embedding/libembedding.a`) and a small demonstration program (`c/embedding/embedding_demo`).

```bash
make
c/embedding/embedding_demo
```

The demo prints the original sequence, the mutated sequence with the payload embedded, and a per-position summary of the encoded bits. Use `make clean` to remove build artefacts.
