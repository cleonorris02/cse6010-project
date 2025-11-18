## Error Detection/Correction Module

A C implementation of the block sum parity scheme described in Na et al. (2020) is available in `c/error_detection`. The module constructs a two-dimensional DNA block augmented with parity nucleotides for each row and column, enabling detection and correction of single-nucleotide mutations in either the data or parity regions.

### Building and running the demo

```
cd c/error_detection
make
./error_detection
```

The demonstration program encrypts three DNA words into a parity-protected block, introduces sample mutations (data, row parity, and column parity), and invokes the recovery routine to restore the correct nucleotides.
